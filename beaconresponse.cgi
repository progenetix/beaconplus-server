#!/usr/bin/perl

# Progenetix & arrayMap site scripts
# Beacon+ server implementation based on arraymap.org | progenetix.org data
# © 2000-2018 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);

use File::Basename;
use JSON;
use List::MoreUtils qw(any apply);
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;
use UUID::Tiny;
use YAML::XS qw(LoadFile);

my $here_path   =   File::Basename::dirname( eval { ( caller() )[1] } );
our $config     =   LoadFile($here_path.'/rsrc/config.yaml') or die print 'Content-type: text'."\n\n¡No config.yaml file in this path!";

=pod

Please see the associated beaconresponse.md

https://beacon.progenetix.test/beaconplus-server/beaconresponse.cgi?datasetId=arraymap&referenceName=9&assemblyId=GRCh38&variantType=DEL&startMin=19,500,000&startMax=21,975,098&endMin=21,967,753&endMax=24,500,000&referenceBases=N&biosamples.biocharacteristics.type.id=ncit:C3224

=cut

#print 'Content-type: text'."\n\n";

if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

# my @names = param();
# foreach (@names) {
#   print Dumper(param($_));
# }
# print $ENV{QUERY_STRING}."\n\n";
# exit;

my $args        =   {};
bless $args;

$args->{datasetPar} =   _getDatasetParams();

# GA4GH variant attributes
$args->{procPar}    =   _getProcessingParams();
$args->{varQpar}    =   _getVariantParams();
$args->{varQpar}    =   _normVariantParams($args->{varQpar});
$args->{varQ}       =   _createVariantQuery($args->{varQpar});
$args->{biosQpar}   =   _getBiosampleParams();
$args->{biosQ}      =   _createBiosampleQuery($args->{biosQpar});

# catching some input errors ###################################################
# TODO: expand ...
$args->{varErrM}    =   _checkVarParams($args->{varQpar});
$args->{biosErrM}   =   _checkBiosParams($args->{biosQpar});
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
  alleleRequest =>  _makePrettyQuery(),
  apiVersion    =>  $config->{api_version},
  url           =>  $config->{url},
  datasetAlleleResponses    =>   $datasetR,
  info          =>  {
    queryString =>  $ENV{QUERY_STRING},
    version     =>  'Beacon+ implementation based on the development branch of the ga4gh-beacon project with custom extensions',
  },

};

print JSON::XS->new->pretty( 1 )->allow_blessed->convert_blessed->canonical()->encode($response)."\n";

exit;

################################################################################
################################################################################
################################################################################

################################################################################
# SUBs #########################################################################
################################################################################

sub _getDatasetParams {

  my $qPar      =   {};

  my %defaults  =   (
    varcoll     =>  'variants',
    callsetcoll =>  'callsets',
  );

  $qPar->{samplecoll}   =   param('samplecoll');
  if ($qPar->{samplecoll} !~ /^\w{3,64}$/) { $qPar->{samplecoll} = 'biosamples' }

  $qPar->{dataset_id}   =   param('datasetId');
  $qPar->{assembly_id}  =   param('assemblyId');
  if ($qPar->{assembly_id} !~ /^\w{3,64}$/) { $qPar->{assembly_id} = 'GRCh38' }

  foreach (keys %defaults) {

    $qPar->{$_} =   param($_);
    if ($qPar->{$_} !~ /^\w{3,64}$/) { $qPar->{$_} = $defaults{$_} }
#    $qPar->{$_} .=   '_'.lc($qPar->{assembly_id});

  }

  $qPar->{dataset_ids}  =   [ param('dataset_ids')];
  if ($qPar->{dataset_ids}->[0] !~ /\w\w\w/) {
    if ($qPar->{dataset_id} =~ /^\w{3,64}$/) { push(@{$qPar->{dataset_ids}}, $qPar->{dataset_id}) }
    else { $qPar->{dataset_ids} = $config->{ dataset_names } }
  }

  return $qPar;

}

################################################################################

sub _getProcessingParams {

  my $qPar      =   {};

  $qPar->{varinfobios}  =   param('phenotypes');
  if ($qPar->{varinfobios} =~ /^(?:y(?:es)?)|1$/i) { $qPar->{varinfobios} = 1 }

  $qPar->{phenotypes}   =   param('phenotypes');
  if ($qPar->{phenotypes} =~ /^(?:y(?:es)?)|1$/i) { $qPar->{phenotypes} = 1 }

  $qPar->{ontologies}   =   [ param('ontologies') ];
  if (
    $qPar->{ontologies}->[0] !~ /^[\w\:]+?$/i
    &&
    $qPar->{phenotypes} > 0
  ) { push(@{$qPar->{ontologies}}, 'NCIT') }
  $qPar->{ontologies}   =   [ grep{ /^[\w\:\/\-]+?$/ } @{$qPar->{ontologies}} ];

  return $qPar;
}

################################################################################

sub _getBiosampleParams {

=pod

Atributes not used (yet):

=cut

  my $qPar      =   {};

#  foreach (qw(
#    id
#    biosamples.biocharacteristics.type.id
#  )) { $qPar->{$_}      =   [ param('biosamples.'.$_) ] }

#  $qPar->{ "id" } = [ param('biosamples.biocharacteristics.type.id') ];
  $qPar->{ "biocharacteristics.type.id" } = [ param('biosamples.biocharacteristics.type.id') ];
#print Dumper($qPar);
  return $qPar;

}

################################################################################

sub _makePrettyQuery {

  my $qPar      =   {};
  foreach my $parameter (qw(
    datasetId
    assemblyId
    referenceName
    referenceBases
    alternateBases
    variantType
    start
    end
    startMin
    startMax
    endMin
    endMax
  )) {
    if (param($parameter) =~ /\w/) {
      $qPar->{$parameter} =   param($parameter);
    }
  };

  # cleaning - e.g. <DEL> => DEL
  foreach my $key (keys %$qPar) {
    $qPar->{$key}   =~ s/[^\w\-\<\>]//g }

  # making sure these are represented as numbers
  foreach my $parameter (qw(
    start
    end
    startMin
    startMax
    endMin
    endMax
  )) {
    if ($qPar->{$parameter} =~ /^\d+?$/) {
      $qPar->{$parameter} *=  1 }
  };

  # TODO: Implement alternate_bases as list

  $qPar->{"biosamples.biocharacteristics.type.id"} = param('biosamples.biocharacteristics.type.id');

  return $qPar;

}

################################################################################

sub _getVariantParams {

=pod

Atributes not used (yet):
  variantset_id
  svlen
  filters_applied
  filters_passed

=cut

  # TODO: Implement alternate_bases as list
  # TODO: Test values based on reference schema
  my $qPar      =   {
    reference_name  =>  param('referenceName') =~ /\w/ ? param('referenceName') : q{},
    reference_bases =>  param('referenceBases') =~ /\w/ ? param('referenceBases') : q{},
    alternate_bases =>  param('alternateBases') =~ /\w/ ? param('alternateBases') : q{},
    variant_type    =>  param('variantType') =~ /\w/ ? param('variantType') : q{},
    start           =>  param('start') =~ /^\d+?$/ ? param('start') : q{},
    end             =>  param('end') =~ /^\d+?$/ ? param('end') : q{},
    start_min       =>  param('startMin') =~ /^[\d\,\']+?$/ ? param('startMin') : q{},
    start_max       =>  param('startMax') =~ /^[\d\,\']+?$/ ? param('startMax') : q{},
    end_min         =>  param('endMin') =~ /[\d\,\']+?$/ ? param('endMin') : q{},
    end_max         =>  param('endMax') =~ /^[\d\,\']+?$/ ? param('endMax') : q{},
  };

  foreach (qw(
    id
  )) { $qPar->{$_}  =   param('variants.'.$_) }

  foreach my $key (keys %$qPar) {
    $qPar->{$key}   =~ s/[^\w\-\.]//g }

  return $qPar;

}

################################################################################

sub _normVariantParams {

  my $qPar      =   $_[0];

  # creating the intervals for range queries, while checking for right order
  # this also fills in min = max if only one parameter has been provided
  # for start or end, respectively
  foreach my $side (qw(start end)) {
    my $parKeys =   [ grep{ /^$side(?:_m(?:(?:in)|(?:ax)))?$/ } keys %$qPar ];
    my $parVals =   [ grep{ /\d/ } @{ $qPar }{ @$parKeys } ];
    $parVals    =   [ apply { $_ =~ s/[\'\,]//g } @{ $parVals } ]; # clean out of formatting characters
    for (@$parVals) { s/[\'\,]//g }
    $parVals    =   [ sort {$a <=> $b} @{ $parVals } ];
    $qPar->{$side.'_range'} =  [ $parVals->[0], $parVals->[-1] ];
  }

  $qPar->{reference_name}   =~  s/chr?o?//i;

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

sub _createBiosampleQuery {

  my $qPar      =   $_[0];
  my @qList;

  foreach my $qKey (keys %{$qPar}) {

    my @thisQlist;

=pod

The query elements are treated through regular expressions, with
    - values being anchored and terminated by either start or end or ':'; this will force complete term matches while disregarding additional prefixes in the data (e.g. "icdom:85003" will match also "pgx:icdom:85003")
    - case insensitive matches

Queries with multiple options for the same attribute are treated as logical "OR".

=cut

    foreach (grep{/.../} @{$qPar->{$qKey}}) { push(@thisQlist, { $qKey => qr/(?:^|\:)$_(?:$|\:)/i }) }

    if (@thisQlist == 1)    { push(@qList, $thisQlist[0]) }
    elsif (@thisQlist > 1)  { push(@qList, {'$or' => [ @thisQlist ] } ) }

  }

=pod

The construction of the query object depends on the detected parameters:

* if empty list => no change, empty object
* if 1 parameter => direct use
* if several parameters are queried => connection through the MongoDB  "$and" constructor

=cut

  if (@qList == 1)    { return $qList[0] }
  elsif (@qList > 1)  { return { '$and' => \@qList } }
  else                { return {} }

}

################################################################################

sub _createVariantQuery {

  my $qPar      =   $_[0];
  my $qObj      =   {};

  #structural query
  if ($qPar->{variant_type} =~ /^D(?:UP)|(?:EL)$/i) {
    $qObj       =   {
      '$and'    => [
        { reference_name      =>  $qPar->{reference_name} },
        { variant_type        =>  $qPar->{variant_type} },
        { start =>  { '$gte'  =>  1 * $qPar->{start_range}->[0] } },
        { start =>  { '$lte'  =>  1 * $qPar->{start_range}->[1] } },
        { end   =>  { '$gte'  =>  1 * $qPar->{end_range}->[0] } },
        { end   =>  { '$lte'  =>  1 * $qPar->{end_range}->[1] } },
      ],
    };

  }

  elsif ($qPar->{variant_type} =~ /^BND$/i) {

    $qObj       =   {
      '$and'    => [
        { reference_name  =>  $qPar->{reference_name} },
        { '$or' =>  [
          { variant_type  =>  'DUP' },
          { variant_type  =>  'DEL' },
        ] },
        { '$or' =>  [
          { '$and'  => [
              { start =>  { '$gte'  =>  1 * $qPar->{start_range}->[0] } },
              { start =>  { '$lte'  =>  1 * $qPar->{start_range}->[1] } },
            ]
          },
          { '$and'  => [
              { end =>  { '$gte'  =>  1 * $qPar->{start_range}->[0] } },
              { end =>  { '$lte'  =>  1 * $qPar->{start_range}->[1] } },
            ]
          },
        ] },
      ],
    };

  }

  # allele query
  elsif ($qPar->{alternate_bases} =~ /^[ATGCN]+?$/) {

    my @qList   =   (
      { reference_name  =>  $qPar->{reference_name} },
      { alternate_bases =>  $qPar->{alternate_bases} },
      { start =>  1 * $qPar->{start} },
    );

    if ($qPar->{reference_bases} =~ /^[ATCG]+?$/) {
      push(
        @qList,
        { reference_bases =>  $qPar->{reference_bases} },
      );
    }

    $qObj       =   { '$and' => \@qList };

  }

  return $qObj;

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

  $MongoDB::Cursor::timeout = 120000;

  my $args      =   shift;
  my $dataset   =   shift;
  my $counts    =   {};
  my $dbCall    =   {};         # recyclable
  my $db        =   $dataset; #.'_ga4gh';
  $db           =~  s/_ga4gh_ga4gh/_ga4gh/;
  $db           =~  s/arraymap_ga4gh/arraymap/;
  my  $dbconn   =   MongoDB::MongoClient->new()->get_database( $db );

=pod

  The ids of biosamples matching (designated) metadata criteria are retrieved. This can be, as in the first example, biosamples with an "characteristic" containing a specific ontology term.

// Dataset's response to a query for information about a specific allele.
message BeaconDatasetAlleleResponse {
  // Identifier of the dataset, as defined in `BeaconDataset`.
  string dataset_id = 1;

  // Indicator of whether the given allele was observed in the dataset.
  //
  // This should be non-null, unless there was an error, in which case
  // `error` has to be non-null.
  bool exists = 2;

  // Dataset-specific error.
  //
  // This should be non-null in exceptional situations only, in which case
  // `exists` has to be null.
  BeaconError error = 3;

  // Frequency of this allele in the dataset. Between 0 and 1, inclusive.
  double frequency = 4;

  // Number of variants matching the allele request in the dataset.
  int64 variant_count = 5;

  // Number of callsets matching the allele request in the dataset.
  int64 call_count = 6;

  // Number of samples matching the allele request in the dataset.
  int64 sample_count = 7;

  // Additional note or description of the response.
  string note = 8;

  // URL to an external system, such as a secured beacon or a system providing
  // more information about a given allele (RFC 3986 format).
  string external_url = 9;

  // A map of additional information.
  map<string, string> info = 13;
}


=cut

  $counts->{variant_count}  =   0;
  $counts->{call_count}     =   0;
  $counts->{frequency}      =   0;
  $counts->{sample_count}   =   0;

  my $bsBioMatchIds     =   []; # biosample ids matching the biosample_Request
  my $csVarBioMatchIds  =   []; # callset ids with varQ and biosQ match
  my $bsVarBioMatchIds  =   []; # biosample ids with varQ and biosQ match

  my $biosBaseNo;
  my $callBaseNo        =   0;
  my $varBaseQ          =   {};
  my $payload;
  my $handover          =   [];

  my $cursor;
  my $vars;

  my $biosAllNo =   $dbconn->get_collection( $args->{datasetPar}->{samplecoll} )->find()->count();

  if (grep{ /.../ } keys %{ $args->{biosQ} } ) {

    $bsBioMatchIds  =   _get_mongo_distinct(
                          $dbconn,
                          $args->{datasetPar}->{samplecoll},
                          'id',
                          $args->{biosQ},
                        );
    $biosBaseNo     =   scalar @$bsBioMatchIds;

    $args->{varQ}   =   { '$and' =>
      [
        $args->{varQ},
        { biosample_id => { '$in' =>  $bsBioMatchIds } },
      ],
    };

  }
  else {
    $biosBaseNo =   $biosAllNo }

  ##############################################################################

  my $fields    =   {
    digest          => 1,
    callset_id      => 1,
    biosample_id    => 1,
    digest          => 1,
  };
  if ($dataset =~ /dgv/i) {
    $fields->{info} =   1 }

  # retrieving all matching variants
  $cursor       =    $dbconn->get_collection( $args->{datasetPar}->{varcoll} )->find( $args->{varQ} )->fields( $fields );
  $vars         =    [ $cursor->all ];

  my %csVarMatches  =   map{ $_->{callset_id} => 1 } @$vars;
  my %bsVarMatches  =   map{ $_->{biosample_id} => 1 } @$vars;
  my %distVars      =   map{ $_->{digest} => 1 } @$vars;

  $counts->{variant_count}  =   scalar keys %distVars;
  if ($counts->{variant_count} > 0) { $counts->{'exists'} = \1 }
#print Dumper( $counts->{variant_count});
#exit;

  my %callSets  =   ();
  $csVarBioMatchIds =   [ keys %csVarMatches ];
  $bsVarBioMatchIds =   [ keys %bsVarMatches ];

  # DGV as a special case ######################################################

  if ($dataset =~ /dgv/i) {

    my %dgsvar  =   ();
    my %dgssvar =   ();

    # In the dgv dataset, "callsets" are in  reality "studies", with the
    # information about samples coming from "sample_size" (of the study)
    # and "count".
    foreach my $var (@$vars) {
      $callSets{ $var->{callset_id} }  =   $var->{info}->{sample_size};
      $counts->{call_count} +=  $var->{info}->{count};
      foreach my $ssvar (split(',', $var->{info}->{supporting_variants})) {
        $dgssvar{$ssvar} =  1;
      }
      foreach my $svar (split(',', $var->{info}->{accession})) {
        $dgsvar{$svar}  =   1;
      }
    }

    $payload    =   {
      supporting_variants   =>   [ sort keys %dgssvar ],
      variants              =>   [ sort keys %dgsvar ],
    };

  }

  # / DGV ######################################################################

  else {

    $counts->{call_count}   =   scalar @$csVarBioMatchIds;
    $counts->{sample_count} =   scalar @$bsVarBioMatchIds;

    if ($biosBaseNo > 0) {
      $counts->{frequency}  =   sprintf "%.4f",  $counts->{sample_count} / $biosBaseNo }

  }

  ##############################################################################

  # storing callset ids for retrieval
  $args->{access_id}    =   create_UUID_as_string();
  my $stored_cs =   {
    _id                 =>  $args->{access_id},
    query_key           =>  'id',
    query_db            =>  $db,
    query_coll          =>  $args->{datasetPar}->{callsetcoll},
    query_values        =>  $csVarBioMatchIds,
  };

  MongoDB::MongoClient->new()->get_database( 'progenetix' )->get_collection( 'querybuffer' )->insert($stored_cs);
  
  push(
    @$handover,
    {
      action    =>  'create CNV histogram from matched callsets',
      url       =>  '/beaconplus-server/beacondeliver.cgi?do=histogram&accessid='.$args->{access_id},
      label     =>  'Histogram'
    },
    {
      action    =>  'export all biosample data of matched callsets',
      url       =>  $config->{url_base}.'/beaconplus-server/beacondeliver.cgi?do=biosamples&accessid='.$args->{access_id},
      label     =>  'Biosamples'
    },
    {
      action    =>  'export all variants of matched callsets',
      url       =>  '/beaconplus-server/beacondeliver.cgi?do=variants&accessid='.$args->{access_id},
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
    query_coll          =>  $args->{datasetPar}->{varcoll},
    query_values        =>  [ keys %distVars ],
  };

  MongoDB::MongoClient->new()->get_database( 'progenetix' )->get_collection( 'querybuffer' )->insert($stored_vars);

  push(
    @$handover,
    {
      action    =>  'retrieve matching variants',
      url       =>  '/beaconplus-server/beacondeliver.cgi?do=variants&accessid='.$args->{varaccess_id},
      label     =>  'Variants'
    }
  );


  ##############################################################################

  if ($dataset =~ /dgv/i) {
    $counts->{frequency}    =   "NA";
    $counts->{sample_count} =   "NA";
  }

  ##############################################################################

  my $bsOntologyTermIds =   _get_mongo_distinct(
                              $dbconn,
                              $args->{datasetPar}->{samplecoll},
                              'biosamples.biocharacteristics.type.id',
                              { id =>  { '$in' => $bsVarBioMatchIds } },
                            );

  ##############################################################################

  my $bsPhenotypeResponse   =   [];

  if ($args->{procPar}->{phenotypes} > 0) {

    foreach my $ontoTerm (@$bsOntologyTermIds) {

      if (grep { $ontoTerm =~ /^$_/i } @{$args->{procPar}->{ontologies}} ) {

        my $ontoNo  =   _count_mongo_distinct(
                          $dbconn,
                          $args->{datasetPar}->{samplecoll},
                          'id',
                          { 'biosamples.biocharacteristics.type.id' => $ontoTerm },
                        );
        my $ontoObs =   _count_mongo_distinct(
                          $dbconn,
                          $args->{datasetPar}->{samplecoll},
                          'id',
                          { '$and' => [
                            { 'biosamples.biocharacteristics.type.id' => $ontoTerm },
                            { id =>  { '$in' => $bsBioMatchIds } },
                          ] },
                        );

        push(
          @$bsPhenotypeResponse,
          {
            id       =>   $ontoTerm,
            count         =>   $ontoNo,
            observations  =>   $ontoObs,
          }
        );

  }}}

  ##############################################################################

  return   {

    datasetId   =>  $dataset,
    "exists"    =>  $counts->{"exists"},
    error       =>  $args->{varErrM}.$args->{biosErrM},
    frequency   =>  $counts->{frequency} * 1,
    variantCount    =>  $counts->{variant_count} * 1,
    callCount   =>  $counts->{call_count} * 1,
    sampleCount =>  $counts->{sample_count} * 1,
    note        =>  ($dataset =~ /dgv/i ? 'Callsets represent the study count.' : q{}),
    externalUrl =>  $config->{url},
    handover    =>  $handover,
    info        =>  {
      callset_access_handle =>  $args->{access_id},
      payload               =>  $payload,
      ontology_selection    =>  $args->{procPar}->{ontologies},
      phenotype_response    =>  $bsPhenotypeResponse,
      description           =>  'The query was against database "'.$db.'", variant collection "'.$args->{datasetPar}->{varcoll}.'". '.$counts->{call_count}.' matched callsets for '.$counts->{variant_count}.' distinct variants. Out of '.$biosAllNo.' biosamples in the database, '.$biosBaseNo.' matched the biosample query; of those, '.$counts->{sample_count}.' had the variant.',
    },

  };

}

################################################################################
# little helpers ###############################################################
################################################################################

sub _get_mongo_distinct {

  my ($dbconn, $collname, $key, $query) =   @_;
  my $dbCall    =   $dbconn->run_command([
                      "distinct"=>  $collname,
                      "key"     =>  $key,
                      "query"   =>  $query,
                    ]);
  return $dbCall->{values};

}

sub _count_mongo_distinct {

  my ($dbconn, $collname, $key, $query) =   @_;
  my $dbCall    =   $dbconn->run_command([
                      "distinct"=>  $collname,
                      "key"     =>  $key,
                      "query"   =>  $query,
                    ]);
  return scalar(@{ $dbCall->{values} } );

}

1;
