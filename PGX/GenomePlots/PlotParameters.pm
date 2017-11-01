package PGX::GenomePlots::PlotParameters;

use Data::Dumper;
use YAML::XS 'LoadFile';

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(read_plot_defaults args_modify_plot_parameters hex2rgb);

################################################################################

sub read_plot_defaults {

  my $plotPars  =   LoadFile('./PGX/rsrc/config/plotdefaults.yaml');  
  return  $plotPars;
  
}

################################################################################

sub args_modify_plot_parameters {

=pod

Expects:
  - the dash prefixed input args
  - the plot parameter object

Returns:
  - modified plot parameter object

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my ($plotPars, $args) =   @_;
  
  # the -plottype | -colorschema specific mappings are processed first & removed 
  # thereafter (there are still fallbacks if no parameters given)
  if (grep{ $args->{'-colorschema'} eq $_ } keys %{ $plotPars->{colorschemas} }) {
    my $colorschema     =   $args->{'-colorschema'};
    foreach (keys %{ $plotPars->{colorschemas}->{ $colorschema } }) {
      if ($plotPars->{colorschemas}->{ $colorschema }->{$_} =~ /^(\#\w{6})$/) {
        $plotPars->{$_} =   $1 }
    }
    delete $plotPars->{colorschemas};
    delete $args->{'-colorschema'};
  }
  if (grep{ $args->{'-plottype'} eq $_ } keys %{ $plotPars->{plottype_values} }) {
    foreach (keys %{ $plotPars->{plottype_values}->{ $args->{'-plottype'} } }) {
      $plotPars->{$_} =   $plotPars->{plottype_values}->{ $args->{'-plottype'} }->{$_} }
    delete $plotPars->{plottype_values};
    delete $args->{'-plottype'};
  }

  foreach my $par (keys %$plotPars) {   
    if ($args->{'-'.$par} =~ /^\-?\#?\w[\.\-\,\w]*?$/) {
      
      # list style parameters are provided comma concatenated
      if (grep{ $par eq $_ } qw(chr2plot label_y_m)) {
        $plotPars->{$par}   =   [ split(',', $args->{'-'.$par}) ] }
      else {
        $plotPars->{$par}   =   $args->{'-'.$par} }   
    }
  }
  
  # derived
  $plotPars->{pixyfactor}   =   1 * $plotPars->{size_plotarea_h_px} / (2 * $plotPars->{value_plot_y_max});
  
  return $plotPars;
  
}

################################################################################

sub hex2rgb {

    my ($r, $g, $b)			=		$_[0]	=~	m/^\#?(\w{2})(\w{2})(\w{2})$/;

    return [ CORE::hex($r), CORE::hex($g), CORE::hex($b) ];

}

1;
