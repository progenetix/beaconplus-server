package beaconPlus::DataExporter;

use Data::Dumper;
use PGX::Helpers::UtilityLibs;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  new
  create_handover_exporter
);

sub new {

=pod

=cut

  my $class     =   shift;
  my $prefetch  =   shift;

  my $self      =   {
    config      =>  $prefetch->{config},
    handover_pre    =>  $prefetch->{handover},
    handover    =>  [],
  };

  bless $self, $class;
  return $self;

}

################################################################################

sub create_handover_exporter {

  my $exporter  =   shift;
    
  foreach (sort keys %{ $exporter->{config}->{handover_types} }) {
  
    my $handoverType  =   $exporter->{config}->{handover_types}->{$_};
    my $handoverPre   =   $exporter->{handover_pre}->{ $handoverType->{handover_method} };
    if ($handoverPre->{target_count} < 1) { next }
  
    push(
      @{ $exporter->{handover} },
      {
        handoverType  =>  {
          id    =>  $handoverType->{id},
          label =>  $handoverType->{label}      
        },
        description =>  $handoverType->{description},
        url     =>  $exporter->{config}->{url_base}.'/beaconplus-server/beacondeliver.cgi?do='.$_.'&accessid='.$handoverPre->{_id},
      }
    );
  
  }

  return $exporter;

}

1;
