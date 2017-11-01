package PGX::GenomePlots::ArrayPlotter;

use Data::Dumper;
use GD::Simple;
use MIME::Base64 qw(encode_base64);
use PGX::GenomePlots::CytobandsPlotter;
use PGX::GenomePlots::PlotParameters;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(return_arrayplot_svg get_arrayplot_area);

################################################################################

sub return_arrayplot_svg {

  our ($plot, $probes, $segments)   =   @_;

  $plot->{Y}            =   $plot->{parameters}->{size_plotmargin_top_px};
  $plot->{areawidth}    =   $plot->{parameters}->{size_plotimage_w_px} - 2 * $plot->{parameters}->{size_plotmargin_px};
  $plot->{basepixfrac}  =   ( $plot->{areawidth} - ($#{ $plot->{parameters}->{chr2plot} } * $plot->{parameters}->{size_chromosome_padding_px}) ) / $plot->{genomesize};
  $plot         =   get_cytobands_svg($plot);
  $plot         =   get_arrayplot_area($plot);
  $plot->{Y}    +=   $plot->{parameters}->{size_plotmargin_bottom_px};
  $plot->{svg}  =   '<svg
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
version="1.1"
id="'.$plot->{plotid}.'"
onload="init(evt)"
width="'.$plot->{parameters}->{size_plotimage_w_px}.'px"
height="'.$plot->{Y}.'px"
style="margin: auto; font-family: Helvetica, sans-serif;">
'.$plot->{svg}.'
</svg>';

  return $plot;

}

sub get_arrayplot_area {

=pod

Expects:
  - the current Y parameter for placing the plot elements on the SVG

Returns:
  - the increased end Y value as start for the next elements

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $areaX_0   =   $plot->{parameters}->{size_plotmargin_px};
  my $areaY_0   =   $plot->{Y} + $plot->{parameters}->{size_plotarea_h_px} / 2;

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", $plot->{referencebounds}->{$refName}->[1] * $plot->{basepixfrac};
		my $areaProbes	=	 [ grep{ $_->{reference_name} eq $refName } @$probes ];
		$areaProbes			=	 [ grep{ $_->{position} <= $plot->{referencebounds}->{$refName}->[1] } @$areaProbes];
		$areaProbes			=	 [ grep{ $_->{position} >= $plot->{referencebounds}->{$refName}->[0] } @$areaProbes ];

    $plot->{svg}    .=  '
<rect x="'.$areaX_0.'" y="'.$plot->{Y}.'" width="'.$areaW.'" height="'.$plot->{parameters}->{size_plotarea_h_px}.'" style="fill: '.$plot->{parameters}->{color_plotarea_hex}.'; fill-opacity: 0.8; " />';

    my $gdDotS      =   1 * $plot->{parameters}->{factor_probedots};
		my $probeArea	  =   GD::Image->new($areaW, $plot->{parameters}->{size_plotarea_h_px}, 1);
    my $gdAreaCol   =   $probeArea->colorAllocate( @{ hex2rgb($plot->{parameters}->{color_plotarea_hex}) } );
    $probeArea->transparent($gdAreaCol);
    my $gdDotcol		=	$probeArea->colorAllocateAlpha(0,0,0,111);
		$probeArea->filledRectangle(0, 0, $areaW, $plot->{parameters}->{size_plotarea_h_px}, $gdAreaCol);

		foreach (@$areaProbes) {
			my $dotX      =   sprintf "%.2f", $plot->{basepixfrac} * $_->{position};
			my $dotY	    =   sprintf "%.2f", ($plot->{parameters}->{size_plotarea_h_px} / 2 - $_->{value} * $plot->{parameters}->{pixyfactor});
      $probeArea->filledEllipse($dotX, $dotY, $gdDotS, $gdDotS, $gdDotcol);
    }

    $plot->{svg}    .=  '
<image
  x="'.$areaX_0.'"
  y="'.$plot->{Y}.'"
  width="'.$areaW.'"
  height="'.$plot->{parameters}->{size_plotarea_h_px}.'"
  xlink:href="data:image/png;base64,'.encode_base64($probeArea->png).'"
/>';

    $areaX_0		+=	$areaW + $plot->{parameters}->{size_chromosome_padding_px};
  }

  $plot->{Y}    +=  $plot->{parameters}->{size_plotarea_h_px};

  return $plot;

}

1;
