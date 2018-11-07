package beaconplus::QueryParameters;

use Data::Dumper;
require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  new
  get_variant_params
);


sub new {

  my $class     =   shift;
  my $args      =   shift;

  my $self      =   {
    variant_param       =>  {}, 
    biosample_param     =>  {}, 
    callset_param       =>  {},
    individual_param    =>  {},
    biosubset_param     =>  {},
    variant_query       =>  {}, 
    biosample_query     =>  {}, 
    callset_query       =>  {},
    individual_query    =>  {},
    biosubset_query     =>  {},
  };

  bless $self, $class;

  return $self;

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

  # TODO: Implement alternate_bases as list
  # TODO: Test values based on reference schema
  my $query     =   shift;
  my $scope     =   'variants';
  
  my $query     =   {
    reference_name  =>  param('referenceName') =~ /\w/ ? param('referenceName') : q{},
    mate_name       =>  param('mateName') =~ /\w/ ? param('mateName') : q{},
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
  )) { $query->{$_}  =   param('variants.'.$_) }

  foreach my $key (keys %$query) {
    $query->{$key}   =~ s/[^\w\-\.]//g }

  return $query;

}