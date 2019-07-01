#!/usr/bin/perl

# Progenetix & arrayMap site scripts
# Beacon+ server implementation based on arraymap.org | progenetix.org data
# © 2000-2019 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);
use Data::Dumper;
use JSON;
use MongoDB;
use MongoDB::MongoClient;
$MongoDB::Cursor::timeout = 36000;

use BeaconPlus::ConfigLoader;
use BeaconPlus::QueryParameters;
use BeaconPlus::QueryExecution;
use BeaconPlus::DataExporter;

=markdown
The "beaconresponse.cgi" script is a server backend for the Beacon protocol,
specifically adapted to connect to _MongoDB_ database backends that adhere to
the core GA4GH model.

#### Example Queries

* <https://beacon.progenetix.org/cgi-bin/beaconresponse.cgi?datasetIds=arraymap&datasetIds=progenetix&referenceName=9&assemblyId=GRCh38&variantType=DEL&startMin=19500000&startMax=21975098&endMin=21967753&endMax=24500000&referenceBases=N&filters=ncit:C3224>
    - focal deletions in the CDKN2A locus, filtered for "ncit:C3224" (Melanoma)
    - this query ia against 2 datasets (progenetix, arraymap)
    - *this query puts some serious strain on the server!*
* <https://beacon.progenetix.org/cgi-bin/beaconresponse.cgi/?datasetIds=arraymap&referenceName=9&assemblyId=GRCh38&variantType=DEL&startMin=19500000&startMax=21975098&endMin=21967753&endMax=24500000&referenceBases=N&filters=icdom-94003>
    - as above, but using "icdom-94003", the internal code for ICD-O 3 "9400/3" (Glioblastoma, NOS)
* <https://beacon.progenetix.org/beaconplus-server/beaconresponse.cgi?datasetIds=dipg&referenceName=17&assemblyId=GRCh38&startMin=7572826&endMax=7579005&referenceBases=*&alternateBases=N&filters=icdot-C71.7&>
    - wildcard range query for allelic variants on chromosome 17, between bases 7572826 and 7579005
    - dataset "dipg"
    - Since the current specification requires a `alternateBases` parameter, and doesn't include wildcards *but* allows "N", the query will __only__ return variants where the replacement has a length of 1. "NN" leads to 2 etc. This may be changed in future versions of the protocol.
* <https://beacon.progenetix.org/beaconplus-server/beaconresponse.cgi?datasetIds=dipg&referenceName=17&assemblyId=GRCh38&start=7577121&referenceBases=G&alternateBases=A&filters=icdot-C71.7&>
		- the "traditional" precise BeaconAlleleRequest, as supported by the first Beacon version
    - dataset "dipg"
=cut

#print 'Content-type: text'."\n\n";

if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

my $config      =   BeaconPlus::ConfigLoader->new();

################################################################################

my $datasetR    =   [];
if (! grep{ /ERROR/i } @{ $config->{query_errors} }) {
  foreach (@{ $config->{ dataset_names } }) {
    push(@$datasetR, _getDataset($config, $_) );
  }
}

my $bcExists    =   \0;
if (grep{ $_->{exists} } @$datasetR) { $bcExists = \1 }
my $response    =   {
  beaconId      =>  $config->{beacon_id},
  exists        =>  $bcExists,
  alleleRequest =>  {
    referenceName	=>	$config->{param}->{referenceName}->[0],
    startMin		=>	1 * $config->{param}->{startMin}->[0],
    start				=>	1 * $config->{param}->{startMin}->[0],
    startMax		=>	1 * $config->{param}->{startMax}->[0],
    endMin			=>	1 * $config->{param}->{endMin}->[0],
    endMax			=>	1 * $config->{param}->{endMax}->[0],
    referenceBases	=>	$config->{param}->{referenceBases}->[0],
    alternateBases	=>	$config->{param}->{alternateBases}->[0],
    variantType			=>	$config->{param}->{variantType}->[0]
  },
  apiVersion    =>  $config->{api_version},
  url           =>  $config->{url},
  datasetAlleleResponses  =>   $datasetR,
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

  my $config    =   shift;
  my $dataset   =   shift;

  # Handover:
  my $handover  =   [];
  my $varResp  	=   [];
  my $varDistNo =   0;

=markdown
#### Beacon Query Execution

Different types of query are run, if parameters for them exist.

The current concept is:

- run queries on the different primary object types (individual, biosample, callset, variant)
- provide the object-level counts as output (and store the ids for handover procedures)
- aggregate on biosamples

* Variant query

* Callset query

* Biosample query

* Individual query

The details of the query execution can be found in the documentation for the
[__BeaconPlus::QueryExecution.pm__](+generated-doc-BeaconPlus-QueryExecution/) module.

=cut

  my $counts    =   {
    frequency   => 0,
    exists      =>  \0,
  };

  my $prefetch  =   BeaconPlus::QueryExecution->new($config, $dataset);
  $prefetch->execute_aggregate_query();

  my $varDistNo  =   $prefetch->{handover}->{'variants::digest'}->{target_count};

  ##############################################################################

  my $exporter  =   BeaconPlus::DataExporter->new($config, $prefetch);
  if ($varDistNo > 0) {
    $exporter->create_handover_exporter();
    $varResp		=		$prefetch->{handover}->{'variants::digest'}->{target_values};
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
    varResp    	=>  $varResp,
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
