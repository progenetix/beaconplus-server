#!/usr/bin/perl

# Progenetix & arrayMap site scripts
# © 2000-2017 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);

use JSON;
use List::Util qw(min max);
use List::MoreUtils qw(any apply);
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;
use UUID::Tiny;

=pod

Please see the associated beaconresponse.md

Query & response processing:

variantQuery (varQ) =>



=cut

if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

my $args        =   {};

$args->{datasetPar}     =   _getDatasetParams();

# GA4GH variant attributes
$args->{procPar}        =   _getProcessingParams();
$args->{varQpar}        =   _getVariantParams();
$args->{varQpar}        =   _normVariantParams($args->{varQpar});
$args->{varQ}           =   _createVariantQuery($args->{varQpar});
$args->{biosQpar}       =   _getBiosampleParams();
$args->{biosQ}          =   _createBiosampleQuery($args->{biosQpar});

# catching some input errors ###################################################

# TODO: expand ...
$args->{errorM}         =   _checkParameters($args->{varQpar});
$args->{queryScope}     =   'datasetAlleleResponses';
$args->{queryType}      =   'alleleRequest';
if ($args->{varQpar}->{variant_type} =~ /^D(?:UP)|(?:EL)$/i) {
  $args->{datasetPar}->{varcoll}        =~  s/_alleles_/_cnv_/;
  $args->{datasetPar}->{callsetcoll}    =~  s/_alleles_/_cnv_/;
}

################################################################################

my $datasetResponses    =   _getDatasetResponses($args);
my $bcExists    =   \0;
if (grep{ $_->{exists} } @$datasetResponses) { $bcExists = \1 }
my $response    =   {
  beacon_id     =>  "org.progenetix:progenetix-beacon",
  exists        =>  $bcExists,
  allele_request        =>  _makePrettyQuery(),
  api_version   =>  "0.4",
  url           =>  'http://progenetix.org/beacon/info/',
  dataset_allele_responses      =>   $datasetResponses,
  info  =>  {
    query_string        =>  $ENV{QUERY_STRING},
    biosample_request   =>  $args->{biosQ},
    version             =>  'Beacon+ implementation based on a development branch of the beacon-team project: https://github.com/ga4gh/beacon-team/pull/94',
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
    varcoll     =>  'variants_alleles',
    callsetcoll =>  'callsets_alleles',
  );

  $qPar->{samplecoll}   =   param('samplecoll');
  if ($qPar->{samplecoll} !~ /^\w{3,64}$/) { $qPar->{samplecoll} = 'biosamples' }

  $qPar->{dataset_id}   =   param('dataset_id');
  $qPar->{assembly_id}  =   param('assembly_id');
  if ($qPar->{assembly_id} !~ /^\w{3,64}$/) { $qPar->{assembly_id} = 'GRCh36' }

  foreach (keys %defaults) {

    $qPar->{$_} =   param($_);
    if ($qPar->{$_} !~ /^\w{3,64}$/) { $qPar->{$_} = $defaults{$_} }
    $qPar->{$_} .=   '_'.lc($qPar->{assembly_id});

  }

  $qPar->{dataset_ids}  =   [ param('dataset_ids')];
  if ($qPar->{dataset_ids}->[0] !~ /\w\w\w/) {
    if ($qPar->{dataset_id} =~ /^\w{3,64}$/) { push(@{$qPar->{dataset_ids}}, $qPar->{dataset_id}) }
    else { $qPar->{dataset_ids} = [ qw(arraymap dipg) ] }
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

  foreach (qw(
    id
    bio_characteristics.ontology_terms.term_id
  )) { $qPar->{$_}      =   [ param('biosamples.'.$_) ] }
#print Dumper($qPar);
  return $qPar;

}

################################################################################

sub _makePrettyQuery {

  my $qPar      =   {};

# TODO: Implement alternate_bases as list

  foreach (qw(
    id
    reference_name
    reference_bases
    alternate_bases
    variant_type
    start
    end
    start_min
    start_max
    end_min
    end_max
  )) { $qPar->{$_}      =   param('variants.'.$_) }

  $qPar->{"bio_characteristics.ontology_terms.term_id"} = param('biosamples.bio_characteristics.ontology_terms.term_id');

  return { '$and' => [ map{ { $_ => $qPar->{$_} } } (grep{ $qPar->{$_} =~ /\w/ } keys %$qPar )] }

}





sub _getVariantParams {

=pod

Atributes not used (yet):
  variant_set_id
  svlen
  filters_applied
  filters_passed

  cipos
  ciend

=cut

  my $qPar      =   {};

# TODO: Implement alternate_bases as list

  foreach (qw(
    id
    reference_name
    reference_bases
    alternate_bases
    variant_type
    start
    end
    start_min
    start_max
    end_min
    end_max
  )) { $qPar->{$_}      =   param('variants.'.$_) }

  foreach (qw(
  )) { $qPar->{$_}      =   [ sort {$a <=> $b } (param('variants.'.$_)) ] }

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
    $parVals    =   [ apply { $_ =~ s/[^\d]//g } @{ $parVals } ];
    $qPar->{$side.'_range'}     =  [ min(@$parVals), max(@$parVals) ];
  }

  $qPar->{reference_name}       =~  s/chr?o?//i;

  return $qPar;

}

################################################################################

sub _checkParameters {

  my $qPar      =   $_[0];
  my $error;

  if ( $qPar->{variant_type} =~ /^D(?:UP)|(?:EL)$/ && ( $qPar->{start_range}->[0] !~ /^\d+?$/ || $qPar->{end_range}->[0] !~ /^\d+?$/ ) ) {
    $error      .=    '"variants.start" (and also start_min, start_max) or "variants.end" (and also end_min, end_max) did not contain a numeric value. ' }
  if ($qPar->{reference_name} !~ /^(?:(?:(?:1|2)?\d)|x|y)$/i) {
    $error      .=    '"variants.reference_name" did not contain a valid value (e.g. "chr17" "8", "X"). ' }
  if ( $qPar->{variant_type} !~ /^D(?:UP)|(?:EL)$/ && $qPar->{alternate_bases} !~ /^[ATGC]+?$/ ) {
    $error      .=    'There was no valid value for either "variants.variant_type" or "variants.alternate_bases". ' }

  return $error;

}

################################################################################

sub _createBiosampleQuery {

  my $qPar      =   $_[0];
  my @qList;

  foreach my $qKey (keys %{$qPar}) {

    my @thisQlist;

    foreach (grep{/.../} @{$qPar->{$qKey}}) { push(@thisQlist, {$qKey => qr/(?:^|\:)$_(?:$|\:)/i}) }

=pod

Queries with multiple options for the same attribute are treated as logical "OR".

=cut

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
  if ($qPar->{variant_type} =~ /^D(?:UP)|(?:EL)$/) {

    $qObj       =   {
      '$and'    => [
        { reference_name        =>  $qPar->{reference_name} },
        { variant_type          =>  $qPar->{variant_type} },
        { start =>  { '$gte'    =>  1 * $qPar->{start_range}->[0] } },
        { start =>  { '$lte'    =>  1 * $qPar->{start_range}->[1] } },
        { end   =>  { '$gte'    =>  1 * $qPar->{end_range}->[0] } },
        { end   =>  { '$lte'    =>  1 * $qPar->{end_range}->[1] } },
      ],
    };

  }

  # allele query
  elsif ($qPar->{alternate_bases} =~ /^[ATGC]+?$/) {

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
  my $db        =   $dataset.'_ga4gh';
  $db           =~  s/_ga4gh_ga4gh/_ga4gh/;
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

  // Number of calls matching the allele request in the dataset.
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
  my $csBioMatchIds     =   []; # callset ids matching the biosample_Request
  my $csVarBioMatchIds  =   []; # callset ids with varQ and biosQ match

  my $biosBaseNo;
  my $callBaseNo        =   0;
  my $varBaseQ          =   {};
  my $payload;

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
    $biosBaseNo =   scalar @$bsBioMatchIds;

    $csBioMatchIds  =   _get_mongo_distinct(
                          $dbconn,
                          $args->{datasetPar}->{callsetcoll},
                          'id',
                          { biosample_id => { '$in' =>  $bsBioMatchIds } },
                        );
    $varBaseQ   =   { 'calls.call_set_id' => { '$in' =>  $csBioMatchIds } };
    $args->{varQ}   =   { '$and' =>
      [
        $args->{varQ},
        $varBaseQ,
      ],
    };

  } else {

    $biosBaseNo =   $biosAllNo;

  }

  ##############################################################################

  # retrieving all matching calls (lists in variants)
  $cursor	      =		$dbconn->get_collection( $args->{datasetPar}->{varcoll} )->find( $args->{varQ} )->fields( { calls => 1 } );
  $vars	        =		[ $cursor->all ];

  $counts->{variant_count}  =   scalar @$vars;
  if ($counts->{variant_count} > 0) { $counts->{'exists'} = \1 }

  my %callSets  =   ();
  my $csIds     =   [];

  # DGV as a special case ######################################################

  if ($dataset =~ /dgv/i) {

    my %dgsvar  =   ();
    my %dgssvar =   ();

    # In the dgv dataset, "callsets" are in  reality "studies", with the
    # information about samples coming from "sample_size" (of the study)
    # and "count".
    foreach my $var (@$vars) {
      foreach (@{ $var->{calls} }) {
        $callSets{ $_->{call_set_id} }  =   $_->{info}->{sample_size};
        $counts->{call_count} +=  $_->{info}->{count};
        foreach my $ssvar (split(',', $_->{info}->{supporting_variants})) {
          $dgssvar{$ssvar}     =   1;
        }
        foreach my $svar (split(',', $_->{info}->{accession})) {
          $dgsvar{$svar}    =   1;
        }
      }
    }

    if ($counts->{c_all} > 0) {
      $counts->{frequency}  =   sprintf "%.4f",  $counts->{call_count} / $counts->{c_all} }

    $payload    =   {
      supporting_variants   =>   [ sort keys %dgssvar ],
      variants              =>   [ sort keys %dgsvar ],
    };

  }

  # / DGV ######################################################################

  else {

    if (grep{ /.../ } keys %{ $args->{biosQ} } ) {
      foreach my $var (@$vars) {
        foreach (@{ $var->{calls} }) {
          my $csid      =   $_->{call_set_id};
          if (grep{ $_ eq $csid } @$csBioMatchIds) {
            $callSets{ $csid }  +=   1;
    }}}}

    else {                              # no csid check
      foreach my $var (@$vars) {
        foreach (@{ $var->{calls} }) {
          $callSets{ $_->{call_set_id} }  +=   1;
    }}}

    $counts->{sample_count} =   _count_mongo_distinct(
                                  $dbconn,
                                  $args->{datasetPar}->{callsetcoll},
                                  'biosample_id',
                                  { id => { '$in' => [ keys %callSets ] } },
                                );

  }

  $csVarBioMatchIds     =   [ keys %callSets ];
  foreach (keys %callSets) { $counts->{call_count} += $callSets{ $_ } }
  if ($biosBaseNo > 0) {
    $counts->{frequency}    =   sprintf "%.4f",  $counts->{sample_count} / $biosBaseNo }


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

  if ($dataset =~ /dgv/i) {
    $counts->{frequency}       =   "NA" }

  ##############################################################################

  my $bsOntologyTermIds =   _get_mongo_distinct(
                              $dbconn,
                              $args->{datasetPar}->{samplecoll},
                              'bio_characteristics.ontology_terms.term_id',
                              { id =>  { '$in' => $bsBioMatchIds } },
                            );

  ##############################################################################

  my $bsPhenotypeResponse       =   [];

  if ($args->{procPar}->{phenotypes} > 0) {

    foreach my $ontoTerm (@$bsOntologyTermIds) {

      if (any { $ontoTerm =~ /^$_/i } @{$args->{procPar}->{ontologies}} ) {

        my $ontoNo  =   _count_mongo_distinct(
                          $dbconn,
                          $args->{datasetPar}->{samplecoll},
                          'id',
                          { 'bio_characteristics.ontology_terms.term_id' => $ontoTerm },
                        );
        my $ontoObs =   _count_mongo_distinct(
                          $dbconn,
                          $args->{datasetPar}->{samplecoll},
                          'id',
                          { '$and' => [
                            { 'bio_characteristics.ontology_terms.term_id' => $ontoTerm },
                            { id =>  { '$in' => $bsBioMatchIds } },
                          ] },
                        );

        push(
          @$bsPhenotypeResponse,
          {
            term_id       =>   $ontoTerm,
            count         =>   $ontoNo,
            observations  =>   $ontoObs,
          }
        );

  }}}

  ##############################################################################

  return   {

    dataset_id          =>  $dataset,
    exists              =>  $counts->{exists},
    error               =>  $args->{errorM},
    frequency           =>  $counts->{frequency},
    variant_count       =>  $counts->{variant_count},
    call_count          =>  $counts->{call_count},
    sample_count        =>  $counts->{sample_count},
    note                =>  ($dataset =~ /dgv/i ? 'Callsets represent the study count.' : q{}),
    external_url        =>  'http://beacon.arraymap.org',
    info                =>  {
      callset_access_handle     =>  $args->{access_id},
      payload                   =>  $payload,
      ontology_selection        =>  $args->{procPar}->{ontologies},
      phenotype_response        =>  $bsPhenotypeResponse,
      description               =>  'The query was against database "'.$db.'", variant collection "'.$args->{datasetPar}->{varcoll}.'". '.$counts->{call_count}.' matched calls for '.$counts->{variant_count}.' distinct variants. Out of '.$biosAllNo.' biosamples in the database, '.$biosBaseNo.' matched the biosample query; of those, '.$counts->{sample_count}.' had the variant.',
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
