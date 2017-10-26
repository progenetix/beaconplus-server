package PGX::GenomePlots::ProbesSegmentsPlotter;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(get_probesplot_svg);

################################################################################

sub get_probesplot_svg {

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

# TODO: evaluate doing this OO, with all standard parameters in the object
  my (
    $plotPars->{CYTOBANDS},
    $plotPars->{GENOINTVS},
    $plotPars->{REFBOUNDS},
    $plotPars->{GENOMESIZE},
    $plotPars,
    $svgId,
    $Y,
    $probes,
    $segments
  )             =   @_;

  my $svg;
  my $plotarea_w  =   $plotPars->{size_plotimage_w_px} - 2 * $plotPars->{size_plotmargin_px};
  my $basePixF  =   ( $plotarea_w - ($#{ $plotPars->{chr2plot} } * $plotPars->{size_chromosome_padding_px}) ) / $plotPars->{GENOMESIZE};
  my $areaX_0   =   $plotPars->{size_plotmargin_px};

  # plot area
  $svg  .=  '
<rect x="'.$plotPars->{size_plotmargin_px}.'" y="'.$Y.'" width="'.$plotarea_w.'" height="'.$plotPars->{size_plotarea_h_px}.'" style="fill: '.$plotPars->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

  $areaX_0      =   $plotPars->{size_plotmargin_px};
  my $areaY_0   =   $Y + $plotPars->{size_plotarea_h_px} / 2;
  foreach my $refName (@{ $plotPars->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", $plotPars->{REFBOUNDS}->{$refName}->[1] * $basePixF;
    my @ind     =  grep{
          $refName eq	$plotPars->{GENOINTVS}->[ $_ ]->{reference_name} } 0..$#{ $plotPars->{GENOINTVS} };







#     foreach my $GL (qw(dupfrequencies delfrequencies)) {
#       $svg	    .=	  '
# <polygon points="'.$areaX_0.' '.$areaY_0;
#       foreach my $i (@ind) {
#         my $X   =		sprintf "%.1f", $areaX_0 + $basePixF * ($plotPars->{GENOINTVS}->[$i]->{end} - ($plotPars->{GENOINTVS}->[$i]->{end} - $plotPars->{GENOINTVS}->[$i]->{start}) / 2);
#         my $H	  =		sprintf "%.1f", $fMaps->{$GL}->[$i] / 2 * $plotPars->{size_plotarea_h_px} / $plotPars->{value_plot_y_max};
#         $svg	  .=	' '.$X.' '.(sprintf "%.1f", ( $GL eq 'delfrequencies' ? $areaY_0 + $H : $areaY_0 - $H) );
#       }
#       $svg	    .=  ' '.(sprintf "%.1f", $areaX_0 + $areaW ).' '.$areaY_0.'" fill="'.($GL =~ /del/i ? $plotPars->{color_var_del_hex} : $plotPars->{color_var_dup_hex}).'" stroke-width="0px" />';
#
#     }




    $areaX_0		+=	$areaW + $plotPars->{size_chromosome_padding_px};
  }

  $Y    +=  $plotPars->{size_plotarea_h_px};

  return ($svg, $Y);

}

1;
