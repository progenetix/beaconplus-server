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

use Data::Dumper;

use beaconPlus::ConfigLoader;
use beaconPlus::QueryParameters;
use beaconPlus::QueryExecution;

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

  $MongoDB::Cursor::timeout = 120000;

  my $dataset   =   shift;
  my $counts    =   {
    variant_count   => 0,
    call_count      => 0,
    sample_count    => 0,
    frequency       => 0,
    exists          =>  \0,
  };

  ##############################################################################

  # Handover:
  my ($prefetch, $buffered);
  my $handover  =   [];

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


  # while not strictly necessary, here the method keys are created as named
  # variables - just for readability
  my $varids_from_aggregation   =   'variants::_id';
  my $csids_from_aggregation    =   'callsets::_id';
  my $biosids_from_aggregation  =   'biosamples::_id';

  $prefetch     =   beaconPlus::QueryExecution->new($config);
  $prefetch->{dataset}  =   $dataset;
  $prefetch->{db_conn}  =   MongoDB::MongoClient->new()->get_database( $dataset );

  my $biosAllNo =   $prefetch->{db_conn}->get_collection( $config->{collection_names}->{biosample_collection} )->find()->count();
  my $biosBaseNo=   $biosAllNo;

  $prefetch->execute_aggregate_query($query);

  ##############################################################################

  if ($ENV{SERVER_NAME} =~ /test/) { $config->{url_base} =  'http://'.$ENV{SERVER_NAME}}
  push(
    @$handover,
    {
      id        =>  'CUSTOM',
      note      =>  'create CNV histogram from matched callsets',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=histogram&accessid='.$prefetch->{handover}->{$csids_from_aggregation}->{_id},
      label     =>  'CNV Histogram'
    },
    {
      id        =>  'CUSTOM',
      note      =>  'export all biosample data of matched callsets',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=biosamples&accessid='.$prefetch->{handover}->{$biosids_from_aggregation}->{_id},
      label     =>  'Biosamples Data'
    },
    {
      id        =>  'CUSTOM',
      note      =>  'export all variants of matched callsets',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=csvariants&accessid='.$prefetch->{handover}->{$csids_from_aggregation}->{_id},
      label     =>  'Callset Variants'
    },
    {
      id        =>  'CUSTOM',
      note      =>  'retrieve matching variants',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=variants&accessid='.$prefetch->{handover}->{$varids_from_aggregation}->{_id},
      label     =>  'Matching Variants'
    }
  );

  ##############################################################################

  $counts->{variant_count}  =   $prefetch->{handover}->{$varids_from_aggregation}->{target_count};
  if ($counts->{variant_count} > 0) { $counts->{'exists'} = \1 }

  $counts->{call_count}   =   $prefetch->{handover}->{$csids_from_aggregation}->{target_count};
  $counts->{sample_count} =   $prefetch->{handover}->{$biosids_from_aggregation}->{target_count};
  if ($biosBaseNo > 0) {
    $counts->{frequency}  =   sprintf "%.4f",  $counts->{sample_count} / $biosBaseNo }

  if (! $counts->{"exists"}) { $handover = [] }

  return  {
    datasetId   =>  $dataset,
    "exists"    =>  $counts->{"exists"},
    error       =>  $query->{query_errors},
    frequency   =>  $counts->{frequency} * 1,
    variantCount    =>  $counts->{variant_count} * 1,
    callCount   =>  $counts->{call_count} * 1,
    sampleCount =>  $counts->{sample_count} * 1,
    note        =>  q{},
    externalUrl =>  $config->{url},
    handover    =>  $handover,
    info        =>  {
      description           =>  'The query was against database "'.$dataset.'", variant collection "'.$config->{collection_names}->{variant_collection}.'". '.$counts->{call_count}.' matched callsets for '.$counts->{variant_count}.' distinct variants. Out of '.$biosAllNo.' biosamples in the database, '.$biosBaseNo.' matched the biosample query; of those, '.$counts->{sample_count}.' had the variant.',
    },

  };

}

1;
