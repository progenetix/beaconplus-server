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
use Data::Dumper;
use UUID::Tiny;
use YAML::XS qw(LoadFile);

my $here_path   =   File::Basename::dirname( eval { ( caller() )[1] } );
our $config     =   LoadFile($here_path.'/rsrc/config.yaml') or die print 'Content-type: text'."\n\n¡No config.yaml file in this path!";

=pod

Please see the associated beaconresponse.md

https://beacon.progenetix.test/beaconplus-server/beaconresponse.cgi?datasetId=arraymap&referenceName=9&assemblyId=GRCh38&variantType=DEL&startMin=19,500,000&startMax=21,975,098&endMin=21,967,753&endMax=24,500,000&referenceBases=N&biosamples.biocharacteristics.type.id=ncit:C3224

https://beacon.progenetix.test/beaconplus-server/beaconresponse.cgi/?datasetId=arraymap&referenceName=9&assemblyId=GRCh38&variantType=DEL&startMin=19,500,000&startMax=21,975,098&endMin=21,967,753&endMax=24,500,000&referenceBases=N&biosamples.biocharacteristics.type.id=icdom:94003
=cut

#print 'Content-type: text'."\n\n";

if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

my $args        =   {};
bless $args;

$args->{datasetPar} =   _getDatasetParams();

use beaconPlus::QueryParameters;
my $query       =   beaconPlus::QueryParameters->new();

# catching some input errors ###################################################
# TODO: expand ...
$args->{varErrM}    =   _checkVarParams($query->{variant_params});
$args->{biosErrM}   =   _checkBiosParams($query->{biosample_params});
$args->{queryScope} =   'datasetAlleleResponses';
$args->{queryType}  =   'alleleRequest';

################################################################################

my $datasetR    =   [];
if ($args->{varErrM} !~ /.../ || $args->{biosErrM} !~ /.../) {
  $datasetR     =    _getDatasetResponses($args) }
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

sub _getDatasetParams {

  my $qPar      =   {};

  $qPar->{dataset_id}   =   param('datasetId');
  $qPar->{assembly_id}  =   param('assemblyId');
  if ($qPar->{assembly_id} !~ /^\w{3,64}$/) { $qPar->{assembly_id} = 'GRCh38' }

  $qPar->{dataset_ids}  =   [ param('dataset_ids')];
  if ($qPar->{dataset_ids}->[0] !~ /\w\w\w/) {
    if ($qPar->{dataset_id} =~ /^\w{3,64}$/) { push(@{$qPar->{dataset_ids}}, $qPar->{dataset_id}) }
    else { $qPar->{dataset_ids} = $config->{ dataset_names } }
  }

  return $qPar;

}

################################################################################

sub _checkVarParams {

  my $qPar      =   $_[0];
  my $errorM;

  if ( $qPar->{variant_type} =~ /^(?:UP)|(?:EL)$/ && ( $qPar->{start_range}->[0] !~ /^\d+?$/ || $qPar->{end_range}->[0] !~ /^\d+?$/ ) ) {
    $errorM     .=    '"startMin" (and also startMax) or "endMin" (and also endMax) did not contain a numeric value - both are required for DUP & DEL. ' }
  if ( $qPar->{variant_type} =~ /^BND$/ && ( $qPar->{start_range}->[0] !~ /^\d+?$/ && $qPar->{end_range}->[0] !~ /^\d+?$/ ) ) {
    $errorM     .=    'Neither "startMin" (and also startMax) or "endMin" (and also endMax) did contain a numeric value - one range is required for BND. ' }
  if ($qPar->{reference_name} !~ /^(?:(?:(?:1|2)?\d)|x|y)$/i) {
    $errorM     .=    '"variants.reference_name" did not contain a valid value (e.g. "chr17" "8", "X"). ' }
  if ( $qPar->{variant_type} !~ /^(?:DUP)|(?:DEL)|(?:BND)$/ && $qPar->{alternate_bases} !~ /^[ATGC]+?$/ ) {
    $errorM     .=    'There was no valid value for either "alternateBases or variantType". ' }

  return $errorM;

}

################################################################################

sub _checkBiosParams {

  my $qPar      =   $_[0];
  my $errorM;

  if ( ! grep{ /\w\w\w/ } @{ $qPar->{"biocharacteristics.type.id"} }) {
    $errorM     .=    'No biosamples.biocharacteristics.type.id was provided.' }

  return $errorM;

}

################################################################################

sub _getDatasetResponses {

  my $args      =   shift;
  my $datasets  =   [];

  foreach (@{ $args->{datasetPar}->{dataset_ids} }) {
    push(@$datasets, _getDataset($args, $_) );
  }

  return $datasets;

}

################################################################################

sub _getDataset {

=pod

=cut

  $MongoDB::Cursor::timeout = 120000;

  my $args      =   shift;
  my $dataset   =   shift;
  my $counts    =   {};
  my $dbCall    =   {}; 
  my $db        =   $dataset;
  my $dbconn    =   MongoDB::MongoClient->new()->get_database( $db );

  $counts->{variant_count}  =   0;
  $counts->{call_count}     =   0;
  $counts->{frequency}      =   0;
  $counts->{sample_count}   =   0;

  my $bsBioMatchIds     =   []; # biosample ids matching the biosample_Request
  my $csVarBioMatchIds  =   []; # callset ids with varant_query and biosample_query match

  my $biosBaseNo;
  my $callBaseNo        =   0;
  my $varBaseQ          =   {};
  my $handover          =   [];

  my $cursor;
  my $vars;

  my $biosAllNo =   $dbconn->get_collection( $config->{collection_names}->{biosample_collection} )->find()->count();

  if (grep{ /.../ } keys %{ $query->{biosample_query} } ) {

    my $distincts   =   $dbconn->run_command([
                          "distinct"=>  $config->{collection_names}->{biosample_collection},
                          "key"     =>  'id',
                          "query"   =>  $query->{biosample_query},
                        ]);
    $bsBioMatchIds  =   $distincts->{values};
    $biosBaseNo     =   @$bsBioMatchIds;

    $query->{variant_query}  =   { '$and' =>
      [
        $query->{variant_query},
        { biosample_id => { '$in' =>  $bsBioMatchIds } },
      ],
    };

  }
  else {
    $biosBaseNo =   $biosAllNo }

  ##############################################################################

  my $varFields =   {
    digest          => 1,
    callset_id      => 1,
    biosample_id    => 1,
  };

  # retrieving all matching variants
  $cursor       =    $dbconn->get_collection( $config->{collection_names}->{variant_collection} )->find( $query->{variant_query} )->fields( $varFields );
  $vars         =    [ $cursor->all ];

  my %csVarMatches  =   map{ $_->{callset_id} => 1 } @$vars;
  my %bsVarMatches  =   map{ $_->{biosample_id} => 1 } @$vars;
  my %distVars      =   map{ $_->{digest} => 1 } @$vars;

  $counts->{variant_count}  =   scalar keys %distVars;
  if ($counts->{variant_count} > 0) { $counts->{'exists'} = \1 }

  $csVarBioMatchIds =   [ keys %csVarMatches ];
  $counts->{call_count}   =   scalar @$csVarBioMatchIds;
  $counts->{sample_count} =   keys %bsVarMatches;

  if ($biosBaseNo > 0) {
    $counts->{frequency}  =   sprintf "%.4f",  $counts->{sample_count} / $biosBaseNo }

  ##############################################################################

  # storing callset ids for retrieval
  $args->{access_id}    =   create_UUID_as_string();
  my $stored_cs =   {
    _id                 =>  $args->{access_id},
    query_key           =>  'id',
    query_db            =>  $db,
    query_coll          =>  $config->{collection_names}->{callset_collection},
    query_values        =>  $csVarBioMatchIds,
  };

  MongoDB::MongoClient->new()->get_database( 'progenetix' )->get_collection( 'querybuffer' )->insert($stored_cs);
  
  push(
    @$handover,
    {
      note      =>  'create CNV histogram from matched callsets',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=histogram&accessid='.$args->{access_id},
      label     =>  'Histogram'
    },
    {
      note      =>  'export all biosample data of matched callsets',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=biosamples&accessid='.$args->{access_id},
      label     =>  'Biosamples'
    },
    {
      note      =>  'export all variants of matched callsets',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=variants&accessid='.$args->{access_id},
      label     =>  'Callsets'
    }
  );

  ##############################################################################

  # storing variant ids for retrieval
  $args->{varaccess_id} =   create_UUID_as_string();
  my $stored_vars  =   {
    _id                 =>  $args->{varaccess_id},
    query_key           =>  'digest',
    query_db            =>  $db,
    query_coll          =>  $config->{collection_names}->{variant_collection},
    query_values        =>  [ keys %distVars ],
  };

  MongoDB::MongoClient->new()->get_database( 'progenetix' )->get_collection( 'querybuffer' )->insert($stored_vars);

  push(
    @$handover,
    {
      note      =>  'retrieve matching variants',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=variants&accessid='.$args->{varaccess_id},
      label     =>  'Variants'
    }
  );

  if (! $counts->{"exists"}) { $handover = [] }
  
  ##############################################################################

  return  {
    datasetId   =>  $dataset,
    "exists"    =>  $counts->{"exists"},
    error       =>  $args->{varErrM}.$args->{biosErrM},
    frequency   =>  $counts->{frequency} * 1,
    variantCount    =>  $counts->{variant_count} * 1,
    callCount   =>  $counts->{call_count} * 1,
    sampleCount =>  $counts->{sample_count} * 1,
    note        =>  q{},
    externalUrl =>  $config->{url},
    handover    =>  $handover,
    info        =>  {
      callset_access_handle =>  $args->{access_id},
      description           =>  'The query was against database "'.$db.'", variant collection "'.$config->{collection_names}->{variant_collection}.'". '.$counts->{call_count}.' matched callsets for '.$counts->{variant_count}.' distinct variants. Out of '.$biosAllNo.' biosamples in the database, '.$biosBaseNo.' matched the biosample query; of those, '.$counts->{sample_count}.' had the variant.',
    },

  };

}

1;
