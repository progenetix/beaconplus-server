package beaconPlus::QueryParameters;

use Data::Dumper;
use CGI qw(param);
require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  new
  read_param_config
  get_variant_params
  norm_variant_params
  create_variant_query
  get_biosample_params
);


sub new {

  use File::Basename;

=pod

    callset_params      =>  {},
    individual_params   =>  {},
    biosubset_params    =>  {},
    callset_query       =>  {},
    individual_query    =>  {},
    biosubset_query     =>  {},

=cut

  my $class     =   shift;
  my $args      =   shift;

  my $self      =   {
    pretty_params   =>  [],
    variant_params  =>  {},
    variant_query   =>  {},
    biosample_params=>  {},
    biosample_query =>  {},
  };

  bless $self, $class;
  
  $self->read_param_config();
  $self->get_variant_params();
  $self->norm_variant_params();
  $self->create_variant_query();
  $self->get_biosample_params();
  $self->create_biosample_query();

  return $self;

}

################################################################################

sub read_param_config {

  use File::Basename;
  use YAML::XS qw(LoadFile);
  
  my $query     =   shift;
  
  my $here_path =   File::Basename::dirname( eval { ( caller() )[1] } ); 
  $query->{param_config}    =   LoadFile($here_path.'/config/query_params.yaml');
  return $query; 

}

################################################################################

sub get_biosample_params {

  my $query     =   shift;
  my $scope     =   'biosamples';

  my $pretty    =   {};

  my $this_par  =   $query->{param_config}->{biosample_params};
  foreach my $q_param (keys %{ $this_par }) {
  
    if ($this_par->{$q_param}->{type} =~ /array/i) {
      if (param($this_par->{$q_param}->{name})) {
        $query->{biosample_params}->{$q_param}  =   [ param($this_par->{$q_param}->{name}) ];
        $pretty->{$this_par->{$q_param}->{name}}=   [ param($this_par->{$q_param}->{name}) ];
      }
    }
    else {
     $query->{biosample_params}->{$q_param}    =   q{};   
      my $val;
      if (param($this_par->{$q_param}->{name})  =~ /./) {
        $val      =   param($this_par->{$q_param}->{name});
        $val      =~  s/[^\w\-\.\:]//g;    
        if ($val  =~ /$this_par->{$q_param}->{pattern}/) {
          $query->{biosample_params}->{$q_param}  =   $val;
          $pretty->{$this_par->{$q_param}->{name}}=   $val;
        }
      }
    }
  }

  push (
    @{ $query->{pretty_params} },
    { $scope    =>  $pretty }
  );
  
  return $query;

}

################################################################################

sub get_variant_params {

=pod

Atributes not used (yet):
  variantset_id
  svlen
  filters_applied
  filters_passed

=cut

  my $query     =   shift;
  my $scope     =   'variants';
  
  my $pretty    =   {};

  my $this_par  =   $query->{param_config}->{variant_params};
  foreach my $q_param (keys %{ $this_par }) {
    $query->{variant_params}->{$q_param}    =   q{};
    my $val;
    if (param($this_par->{$q_param}->{name})  =~ /./) {
      $val      =   param($this_par->{$q_param}->{name});
      $val      =~  s/[^\w\-\.]//g;
    
      if ($val  =~ /$this_par->{$q_param}->{pattern}/) {
        if ($this_par->{$q_param}->{type} =~/(?:num)|(?:int)|(?:float)/i) {
          $pretty->{ $this_par->{$q_param}->{name} }  =   1 * $val;
          $query->{variant_params}->{$q_param}        =   1 * $val;
        }
        else {
          $pretty->{ $this_par->{$q_param}->{name} }  =   $val;
          $query->{variant_params}->{$q_param}        =   $val;
        }
      }
    }
  }
  
  push (
    @{ $query->{pretty_params} },
    { $scope    =>  $pretty }
  );

  return $query;

}

################################################################################

sub norm_variant_params {

  my $query     =   shift;

  # creating the intervals for range queries, while checking for right order
  # this also fills in min = max if only one parameter has been provided
  # for start or end, respectively
  foreach my $side (qw(start end)) {
    my $parKeys =   [ grep{ /^$side(?:_m(?:(?:in)|(?:ax)))?$/ } keys %{ $query->{variant_params} } ];
    my $parVals =   [ grep{ /\d/ } @{ $query->{variant_params} }{ @$parKeys } ];
    $parVals    =   [ sort {$a <=> $b} @{ $parVals } ];
    $query->{variant_params}->{$side.'_range'} =  [ $parVals->[0], $parVals->[-1] ];
  }

  $query->{variant_params}->{reference_name}   =~  s/chr?o?//i;

  return $query;

}

################################################################################

sub create_variant_query {

  my $query     =   shift; 

  #structural query
  if ($query->{variant_params}->{variant_type} =~ /^D(?:UP)|(?:EL)$/i) {
    $query->{variant_query}       =   {
      '$and'    => [
        { reference_name      =>  $query->{variant_params}->{reference_name} },
        { variant_type        =>  $query->{variant_params}->{variant_type} },
        { start =>  { '$gte'  =>  1 * $query->{variant_params}->{start_range}->[0] } },
        { start =>  { '$lte'  =>  1 * $query->{variant_params}->{start_range}->[1] } },
        { end   =>  { '$gte'  =>  1 * $query->{variant_params}->{end_range}->[0] } },
        { end   =>  { '$lte'  =>  1 * $query->{variant_params}->{end_range}->[1] } },
      ],
    };

  }

  elsif ($query->{variant_params}->{variant_type} =~ /^BND$/i) {

    $query->{variant_query}       =   {
      '$and'    => [
        { reference_name  =>  $query->{variant_params}->{reference_name} },
        { '$or' =>  [
          { variant_type  =>  'DUP' },
          { variant_type  =>  'DEL' },
        ] },
        { '$or' =>  [
          { '$and'  => [
              { start =>  { '$gte'  =>  1 * $query->{variant_params}->{start_range}->[0] } },
              { start =>  { '$lte'  =>  1 * $query->{variant_params}->{start_range}->[1] } },
            ]
          },
          { '$and'  => [
              { end =>  { '$gte'  =>  1 * $query->{variant_params}->{start_range}->[0] } },
              { end =>  { '$lte'  =>  1 * $query->{variant_params}->{start_range}->[1] } },
            ]
          },
        ] },
      ],
    };

  }

  # allele query
  elsif ($query->{variant_params}->{alternate_bases} =~ /^[ATGCN]+?$/) {

    my @qList   =   (
      { reference_name  =>  $query->{variant_params}->{reference_name} },
      { alternate_bases =>  $query->{variant_params}->{alternate_bases} },
      { start =>  1 * $query->{variant_params}->{start} },
    );

    if ($query->{variant_params}->{reference_bases} =~ /^[ATCG]+?$/) {
      push(
        @qList,
        { reference_bases =>  $query->{variant_params}->{reference_bases} },
      );
    }

    $query->{variant_query}       =   { '$and' => \@qList };

  }

  return $query;

}

################################################################################

sub create_biosample_query {

  my $query     =   shift;

  my @qList;

  foreach my $qKey (keys %{ $query->{biosample_params} }) {

    my @thisQlist;

=pod
Queries with multiple options for the same attribute are treated as logical "OR".
=cut

    if (ref $query->{biosample_params}->{$qKey} eq 'ARRAY') {
      foreach (@{ $query->{biosample_params}->{$qKey} }) { push(@thisQlist, { $qKey => qr/^$_/i }) }
    } else {
      push(@thisQlist, { $qKey => qr/^$query->{biosample_params}->{$qKey}/i }),
    }

    if (@thisQlist == 1)    { push(@qList, $thisQlist[0]) }
    elsif (@thisQlist > 1)  { push(@qList, {'$or' => [ @thisQlist ] } ) }

  }

=pod

The construction of the query object depends on the detected parameters:

* if empty list => no change, empty object
* if 1 parameter => direct use
* if several parameters are queried => connection through the MongoDB  "$and" constructor

=cut

  if (@qList == 1)    { $query->{biosample_query} =   $qList[0] }
  elsif (@qList > 1)  { $query->{biosample_query} =   { '$and' => \@qList } }

  return $query;
  
}


1;