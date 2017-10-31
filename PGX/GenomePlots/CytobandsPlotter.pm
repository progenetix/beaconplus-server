package PGX::GenomePlots::CytobandsPlotter;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(get_cytobands_svg);

################################################################################

sub get_cytobands_svg {

=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $plot      =   shift;

  my $areaX_0   =   $plot->{parameters}->{size_plotmargin_px};

  $plot->{svg}  .=  '
  <defs>';
  $plot->{svg}  .=  _cytoband_svg_gradients($plot->{plotid});
  $plot->{svg}  .=  '
  <style type="text/css">
  <![CDATA[
  .chrlab { text-anchor: middle; font-size: '.$plot->{parameters}->{size_text_px}.'px }
  ]]>
  </style>';
  $plot->{svg}  .=  '
  </defs>';

  my $chrolabY  =   $plot->{Y} + $plot->{parameters}->{size_text_px} + $plot->{parameters}->{size_chromosome_padding_px};
  my $chroBandY =   $chrolabY + $plot->{parameters}->{size_chromosome_padding_px};

  foreach my $refName (@{ $plot->{parameters}->{chr2plot} }) {

    my $areaW   =  sprintf "%.1f", $plot->{referencebounds}->{$refName}->[1] * $plot->{basepixfrac};
    my $chroX   =  sprintf "%.1f", $areaX_0 + $areaW / 2;

    $plot->{svg}        .=  '
<text x="'.$chroX.'" y="'.$chrolabY.'" class="chrlab">'.$refName.'</text>';

    my $plotInd =  [ grep{
          $refName eq	$plot->{cytobands}->[$_]->{reference_name} } 0..$#{$plot->{cytobands}} ];
    foreach my $i (@$plotInd) {
      my $cbX	  =		sprintf "%.1f", $areaX_0 + $plot->{basepixfrac} * $plot->{cytobands}->[$i]->{start};
      my $cbW	  =		sprintf "%.1f", $plot->{basepixfrac} * ($plot->{cytobands}->[$i]->{end} - $plot->{cytobands}->[$i]->{start});
      my $cbY	  =		$chroBandY;
      my $cbH		=   $plot->{parameters}->{size_chromosome_w_px};

      # terminal and centromere bands are more narrow
      if (
        ($plot->{cytobands}->[$i]->{stain} =~ /cen/i)
        ||
        $plot->{cytobands}->[$i]->{end} >= ($plot->{referencebounds}->{$refName}->[1] - 500)
        ||
        $plot->{cytobands}->[$i]->{start} <= ($plot->{referencebounds}->{$refName}->[0] + 500) ) {
        $cbY	  +=	(0.1*$plot->{parameters}->{size_chromosome_w_px});
        $cbH	  -=	(0.2*$plot->{parameters}->{size_chromosome_w_px});
      }

      # acrocentric "stalk" bands are slim
      if ($plot->{cytobands}->[$i]->{stain} =~ /stalk/i) {
        $cbY    +=	(0.3*$plot->{parameters}->{size_chromosome_w_px});
        $cbH	  -=	(0.6*$plot->{parameters}->{size_chromosome_w_px});
      }

      $plot->{svg}      .=  '
<rect x="'.$cbX.'" y="'.$cbY.'" width="'.$cbW.'" height="'.$cbH.'" style="fill: url(#'.$plot->{plotid}.$plot->{cytobands}->[$i]->{stain}.'); " />';

    }
    $areaX_0	  +=	$areaW + $plot->{parameters}->{size_chromosome_padding_px};
  }

  $plot->{Y}    =  $chroBandY + $plot->{parameters}->{size_chromosome_w_px} + $plot->{parameters}->{size_chromosome_padding_px};

  return $plot;

}

################################################################################

sub _cytoband_svg_gradients {

	my $id			  =	  $_[0];

	return '
<linearGradient id="'.$id.'gpos100" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="rgb(39,39,39)" />
		<stop offset="100%" stop-color="rgb(0,0,0)" />
</linearGradient>
<linearGradient id="'.$id.'gpos75" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="rgb(87,87,87)" />
		<stop offset="100%" stop-color="rgb(39,39,39)" />
</linearGradient>
<linearGradient id="'.$id.'gpos50" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="rgb(196,196,196)" />
		<stop offset="100%" stop-color="rgb(111,111,111)" />
</linearGradient>
<linearGradient id="'.$id.'gpos25" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="rgb(223,223,223)" />
		<stop offset="100%" stop-color="rgb(196,196,196)" />
</linearGradient>
<linearGradient id="'.$id.'gneg" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="white" />
		<stop offset="100%" stop-color="rgb(223,223,223)" />
</linearGradient>
<linearGradient id="'.$id.'gvar" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="rgb(196,196,196)" />
		<stop offset="100%" stop-color="rgb(111,111,111)" />
</linearGradient>
<linearGradient id="'.$id.'stalk" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="rgb(39,39,39)" />
		<stop offset="100%" stop-color="rgb(0,0,0)" />
</linearGradient>
<linearGradient id="'.$id.'acen" x1="0%" x2="0%" y1="0%" y2="80%" spreadMethod="reflect">
		<stop offset="0%" stop-color="rgb(163,55,247)" />
		<stop offset="100%" stop-color="rgb(138,43,226)" />
</linearGradient>
';

}


1;
