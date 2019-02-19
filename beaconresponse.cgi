#!/usr/bin/perl

# Progenetix & arrayMap site scripts
# Beacon+ server implementation based on arraymap.org | progenetix.org data
# © 2000-2018 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);

use File::Basename;
use JSON;
use MongoDB;
use MongoDB::MongoClient;
$MongoDB::Cursor::timeout = 36000;

use Data::Dumper;

use beaconPlus::ConfigLoader;
use beaconPlus::QueryParameters;
use beaconPlus::QueryExecution;
use beaconPlus::DataExporter;

my $config      =   beaconPlus::ConfigLoader->new();
my $datasets    =   [ param('datasetIds') ];
if ($datasets->[0] =~ /\w\w\w/) {
  $config->{ dataset_names }  =   [];
  foreach (grep{ /\w\w\w/ } @$datasets) {
    push(@{ $config->{ dataset_names } }, $_);
  }
}

=pod

Please see the associated beaconresponse.md

https://beacon.progenetix.test/beaconplus-server/beaconresponse.cgi?datasetIds=arraymap&referenceName=9&assemblyId=GRCh38&variantType=DEL&startMin=19,500,000&startMax=21,975,098&endMin=21,967,753&endMax=24,500,000&referenceBases=N&biosamples.biocharacteristics.type.id=ncit:C3224
https://beacon.progenetix.test/beaconplus-server/beaconresponse.cgi/?datasetIds=arraymap&referenceName=9&assemblyId=GRCh38&variantType=DEL&startMin=19,500,000&startMax=21,975,098&endMin=21,967,753&endMax=24,500,000&referenceBases=N&biosamples.biocharacteristics.type.id=icdom:94003

=cut

#print 'Content-type: text'."\n\n";

if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }
 
my $query       =   beaconPlus::QueryParameters->new($config);
#print Dumper($query);
################################################################################

my $datasetR    =   [];
if (! grep{ /ERROR/i } @{ $query->{query_errors} }) {
  foreach (@{ $config->{ dataset_names } }) {
    push(@$datasetR, _getDataset($_) );
  }
}

my $bcExists    =   \0;
if (grep{ $_->{exists} } @$datasetR) { $bcExists = \1 }
my $response    =   {
  beaconId      =>  $config->{beacon_id},
  exists        =>  $bcExists,
  alleleRequest =>  $query->{pretty_params},
  apiVersion    =>  $config->{api_version},
  url           =>  $config->{url},
  datasetAlleleResponses    =>   $datasetR,
  info          =>  {
    queryString =>  $ENV{QUERY_STRING},
    version     =>  'Beacon+ implementation based on the development branch of the ELIXIR Beacon project with custom extensions',
  },
};

print JSON::XS->new->pretty( 1 )->encode($response)."\n";

exit;

################################################################################
################################################################################
# SUBs #########################################################################
################################################################################
################################################################################

sub _getDataset {

=pod

=cut

  my $dataset   =   shift;

 ##############################################################################

  # Handover:
  my $handover  =   [];
  my $variantResponses  = [];
  my $varDistinctCount  =   0;

=pod

Different types of query are run, if parameters for them exist.

The current concept is:

- run queries on the different primary object types (individual, biosample, callset, variant)
- provide the object-level counts as output (and store the ids for handover procedures)
- aggregate on biosamples

### Variant query

### Callset query

This is just an aggregator, since callsets are currently just wrapper objects (experiments delivering sets of variants).

### Biosample query

### Individual query

- e.g. species, sex ...

=cut
  my $counts    =   {
    frequency       => 0,
    exists          =>  \0,
  };

  my $prefetch  =   beaconPlus::QueryExecution->new($config);
  $prefetch->{dataset}  =   $dataset;
  $prefetch->{db_conn}  =   MongoDB::MongoClient->new()->get_database( $prefetch->{dataset} );
  $prefetch->execute_aggregate_query($query);
  
  ##############################################################################

  if ($ENV{SERVER_NAME} =~ /test/) { $config->{url_base} =  'http://'.$ENV{SERVER_NAME}}

  my $exporter  =   beaconPlus::DataExporter->new($prefetch);
  if ($prefetch->{counts}->{'variants_match_count'} > 0) {
    $exporter->create_handover_exporter();
    $prefetch->aggregate_variants();
    $variantResponses =   $prefetch->{variantResponses};
  }

  ##############################################################################

  if ($prefetch->{counts}->{'variants_match_count'} > 0) { $counts->{'exists'} = \1 }

  if ($prefetch->{counts}->{'biosamples_base_count'} > 0) {
    $counts->{frequency}  =   sprintf "%.4f",  $prefetch->{counts}->{'biosamples_match_count'} / $prefetch->{counts}->{'biosamples_base_count'} }

  if (! $counts->{"exists"}) { $handover = [] }

  return  {
    datasetId   =>  $dataset,
    "exists"    =>  $counts->{"exists"},
    error       =>  $prefetch->{error},
    frequency   =>  $counts->{frequency} * 1,
    variantCount    =>  $prefetch->{counts}->{variants_match_count} * 1,
    variantResponses    =>  $variantResponses,
    callCount   =>  $prefetch->{counts}->{callsets_match_count} * 1,
    sampleCount =>  $prefetch->{counts}->{biosamples_match_count} * 1,
    note        =>  q{},
    externalUrl =>  $config->{url},
    datasetHandover =>  $exporter->{handover},
    info        =>  {
      description   =>  'The query was against database "'.$dataset.'", variant collection "'.$config->{collection_names}->{variant_collection}.'". '.$prefetch->{counts}->{callsets_match_count}.' matched callsets for '.$prefetch->{counts}->{variants_match_count}.' distinct variants. Out of '.$prefetch->{counts}->{biosamples_base_count}.' biosamples in the database, '.$prefetch->{counts}->{biosamples_query_match_count}.' matched the biosample query; of those, '.$prefetch->{counts}->{biosamples_match_count}.' had the variant.',
      distinctVarCount  =>  $prefetch->{handover}->{'variants::digest'}->{target_count},
    },

  };

}

1;
