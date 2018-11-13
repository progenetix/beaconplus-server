package beaconPlus::QueryExecution;

use Data::Dumper;
use UUID::Tiny;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  new
  prefetch_data
  create_handover_object
);

sub new {

  use File::Basename;

=pod

=cut

  my $class     =   shift;
  my $config    =   shift;

  my $self      =   {
    config      =>  $config,
    handover_coll   =>   MongoDB::MongoClient->new()->get_database( $config->{handover_db} )->get_collection( $config->{handover_coll} ),
  };

  bless $self, $class;
  return $self;

}

################################################################################

sub prefetch_data {

  my $prefetch  =   shift;
  my $method    =   shift;
  my $query     =   shift;

  $prefetch->{dataset}  ||= 'arraymap';

  my (
    $source_c,
    $source_k,
    $target_c,
    $target_k
  )             =   split('::', $method);
  
  # two components are interpolated to same output
  if ($target_c !~  /\w/) { $target_c = $source_c }
  if ($target_k !~  /\w/) { $target_k = $source_k }
  
  my $distincts =   MongoDB::MongoClient->new()->get_database( $prefetch->{dataset} )->run_command([
                      "distinct"=>  $source_c,
                      "key"     =>  $source_k,
                      "query"   =>  $query,
                    ]);

  my $distVals  =   $distincts->{values};
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

sub create_handover_object {

  my $prefetch  =   shift;
  my $method    =   shift;
  my $query     =   shift;

  $prefetch->prefetch_data($method, $query);    
  $prefetch->{handover}->{$method}->{_id} =   create_UUID_as_string();

  $prefetch->{handover_coll}->insert( $prefetch->{handover}->{$method} );
    
  return $prefetch;

}

1;