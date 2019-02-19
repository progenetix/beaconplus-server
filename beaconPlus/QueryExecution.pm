package beaconPlus::QueryExecution;

use Data::Dumper;
use UUID::Tiny;
use PGX::Helpers::UtilityLibs;
use beaconPlus::ConfigLoader;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  new
  prefetch_data
  create_handover_object
  execute_aggregate_query
);

sub new {

  use File::Basename;
  $MongoDB::Cursor::timeout = 120000;

=pod

=cut

  my $class     =   shift;
  my $config    =   shift;

  my $self      =   {
    config      =>  $config,
    handover_coll   =>   MongoDB::MongoClient->new()->get_database( $config->{handover_db} )->get_collection( $config->{handover_coll} ),
    error       =>  [],
  };

  bless $self, $class;
  return $self;

}

################################################################################

sub prefetch_data {

  my $prefetch  =   shift;
  my $method    =   shift;
  my $query     =   shift;
  my $filters   =   shift;
  
  $prefetch->{dataset}  ||= 'arraymap';
  $prefetch->{filters}  ||= {};
  my (
    $source_c,
    $source_k,
    $target_c,
    $target_k
  )             =   split('::', $method);

  # two components are interpolated to same output
  if ($target_c !~  /\w/) { $target_c = $source_c }
  if ($target_k !~  /\w/) { $target_k = $source_k }

  my $distincts =   $prefetch->{db_conn}->run_command([
                      "distinct"=>  $source_c,
                      "key"     =>  $source_k,
                      "query"   =>  $query,
                    ]);

  my $distVals  =   $distincts->{values};
  if ($prefetch->{filters}->{randno} =~ /^\d+?$/) {
    $distVals   =   RandArr($distVals, $prefetch->{filters}->{randno}) }
  my $distCount =   scalar @{ $distVals };

  $prefetch->{handover}->{$method}  =   {
    source_db         =>  $prefetch->{dataset},
    source_collection =>  $source_c,
    source_key        =>  $source_k,
    target_collection =>  $target_c,
    target_key        =>  $target_k,
    target_values     =>  $distVals,
    target_count      =>  $distCount,
  };

  return $prefetch;

}

################################################################################

sub get_base_counts {

  my $prefetch  =   shift;
  
  $prefetch->{dataset}  ||= 'arraymap';
  $prefetch->{filters}  ||= {};

  $prefetch->{counts}   =   {};
  foreach (qw(callsets biosamples variants individuals)) {
    $prefetch->{counts}->{$_.'_base_count'}  =   $prefetch->{db_conn}->get_collection($_)->find()->count();
  }

  return $prefetch;

}

################################################################################

sub create_handover_object {

  my $prefetch  =   shift;
  my $method    =   shift;
  my $query     =   shift;

  $prefetch->prefetch_data($method, $query);
  $prefetch->{handover}->{$method}->{_id} =   create_UUID_as_string();

  $prefetch->{handover_coll}->insert( $prefetch->{handover}->{$method} );

  return $prefetch;

}

################################################################################

sub execute_aggregate_query {

  my $prefetch  =   shift;
  my $query     =   shift;
  
  $prefetch->get_base_counts();
  push(
    @{$prefetch->{error}},
    @{$query->{query_errors}}  
  );

  # the main prefetch method here is the retrieval of the 
  # "callsets" collection's "id" values
  my $csidHO    =   'callsets::id';
  
  if (! grep{ /.../ }
    keys %{ $query->{callset_query} },
    keys %{ $query->{queries}->{biosamples} },
    keys %{ $query->{queries}->{variants} }
  ) { return }

# 1. Checking for a callsets query & return iq query but no matches
  if (grep{ /.../ } keys %{ $query->{callset_query} } ) {
    $prefetch->prefetch_data( 'callsets::id', $query->{callset_query} );
    if ($prefetch->{handover}->{'callsets::id'}->{target_count} < 1) { 
      return $prefetch } 
  }
  
# 2. Checking for a biosamples query & return if query but no matches
  if (grep{ /.../ } keys %{ $query->{queries}->{biosamples} } ) {
    my $thisM   =   'biosamples::id';
    $prefetch->prefetch_data( $thisM, $query->{queries}->{biosamples} );
    $prefetch->{counts}->{biosamples_query_match_count}  =   $prefetch->{handover}->{$thisM}->{target_count};
    if ($prefetch->{handover}->{$thisM}->{target_count} < 1) { 
      return $prefetch }
# 3. If biosamples matches, retrieve the callset_id values; if there had been 
#    matches before => intersect through added "$in" query
    my $thisQ   =   { 'biosample_id' => { '$in' =>  $prefetch->{handover}->{'biosamples::id'}->{target_values} } };
    if ($prefetch->{handover}->{'callsets::id'}->{target_count} > 0) {
      $thisQ      =   { '$and' =>
        [
          $thisQ,
          { 'id' => { '$in' =>  $prefetch->{handover}->{'callsets::id'}->{target_values} } },
        ],
      };
    }
    $prefetch->prefetch_data('callsets::id', $thisQ);
    if ($prefetch->{handover}->{'callsets::id'}->{target_count} < 1) { 
      return $prefetch } 
  }

# 4. Checking for a variants query; if there had been callset matches before
#    => intersect through added "$in" query; return iq query but no matches
  if (grep{ /.../ } keys %{ $query->{queries}->{variants} } ) {
    my $thisQ   =   $query->{queries}->{variants};
    if ($prefetch->{handover}->{'callsets::id'}->{target_count} > 0) {

      $thisQ    =   { '$and' =>
        [
          $thisQ,
           { 'callset_id' => { '$in' =>  $prefetch->{handover}->{'callsets::id'}->{target_values} } },
        ],
      };
    }
    $prefetch->create_handover_object('variants::_id', $thisQ);

    if ($prefetch->{handover}->{'variants::_id'}->{target_count} < 1) { 
      return $prefetch }

    # after all variants have been identified, their _id values are used
    # to retrieve other associated data (digest, callset_id)
    $thisQ      =    { '_id' => { '$in' => $prefetch->{handover}->{'variants::_id'}->{target_values} } };
    $prefetch->prefetch_data('variants::callset_id', $thisQ);
    
    # just to overwrite the key with the standard $method
    $prefetch->{handover}->{'callsets::id'}->{target_values} =   $prefetch->{handover}->{'variants::callset_id'}->{target_values};
    $prefetch->{handover}->{'callsets::id'}->{target_count}  =   $prefetch->{handover}->{'variants::callset_id'}->{target_count};

    $prefetch->prefetch_data('variants::digest', $thisQ);

  }

# Up to here, queries against callsets, biosamples and variants have all been reduced to the callsets::id values.

# 6. Reset the biosample handover object & store it, from the current callsets

  $prefetch->prefetch_data(
    'callsets::biosample_id',
    { 'id' => { '$in' => $prefetch->{handover}->{'callsets::id'}->{target_values} } },
  );  
  $prefetch->create_handover_object(
    'biosamples::_id',
    { 'id' => { '$in' => $prefetch->{handover}->{'callsets::biosample_id'}->{target_values} } },
  );
  $prefetch->prefetch_data(
    'biosamples::individual_id',
    { '_id' => { '$in' => $prefetch->{handover}->{'biosamples::_id'}->{target_values} } },
  ); 
  $prefetch->create_handover_object(
    'individuals::_id',
    { 'id' => { '$in' => $prefetch->{handover}->{'biosamples::individual_id'}->{target_values} } },
  );
  $prefetch->create_handover_object(
    'callsets::_id',
    { 'id' => { '$in' => $prefetch->{handover}->{'callsets::id'}->{target_values} } },
  );
  
  foreach (qw(callsets biosamples variants individuals)) {
   if ($prefetch->{handover}->{$_.'::_id'}->{target_count}) {
      $prefetch->{counts}->{$_.'_match_count'} =   $prefetch->{handover}->{$_.'::_id'}->{target_count}; 
    }
  }
  if ($prefetch->{handover}->{'variants::digest'}->{target_count}) {
    $prefetch->{counts}->{'variants_distinct_match_count'} =   $prefetch->{handover}->{$_.'::_id'}->{target_count}; 
  }

  return $prefetch;

}

################################################################################

sub aggregate_variants {

  my $prefetch  =   shift;
  $prefetch->{variantResponses} =   [];
  
  if ($prefetch->{handover}->{'variants::digest'}->{target_count} < 1) { 
    return $prefetch }
        
  if ($prefetch->{handover}->{'variants::digest'}->{target_count} > $prefetch->{config}->{max_distinct_variants}) { 
    $prefetch->{variantResponses} =   $prefetch->{handover}->{'variants::digest'}->{target_values};
    push(
      @{$prefetch->{error}},
      'WARNING: More than '.$prefetch->{config}->{max_distinct_variants}.' distict variants => no variant statistics will be performed and only digests are being listed in "variantResponses".'
    );
    return $prefetch;
  }
  
  foreach (@{ $prefetch->{handover}->{'variants::digest'}->{target_values} }) {
    
    my $var     =   $prefetch->{db_conn}->get_collection('variants')->find_one( { 'digest' => $_ } );
    my $parsed  =   {
      referenceName   =>  $var->{reference_name},
      start     =>  $var->{start}->[0],
      referenceBases  =>  $var->{reference_bases},
      alternateBases  =>  $var->{alternate_bases}->[0],
      count     =>  $prefetch->{db_conn}->get_collection('variants')->find({ 'digest' => $_ })->count(),
    };
    
    if ($var->{end}->[-1] > $var->{start}->[0]) {
      $parsed->{end}    =   $var->{end}->[-1] }
    if ($var->{variant_type} =~ /.../) {
      $parsed->{variantType}  =   $var->{variant_type} }
    
    push(
      @{ $prefetch->{variantResponses} },
      $parsed
    );
  
  }

  return $prefetch;

}

1;
