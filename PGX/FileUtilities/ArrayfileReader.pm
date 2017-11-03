package PGX::FileUtilities::ArrayfileReader;

use Data::Dumper;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(read_probefile read_segmentfile);

################################################################################

sub read_probefile {

=pod

Expects:
  - a standard Progenetix style probe file

	ID	chro	pos	log2
	cnvi0111187	17	35295593	0.0859121900
	cnvi0111188	8	65499402	-0.1438023000
	cnvi0111189	2	177061178	-0.0113166000
	cnvi0111190	5	70255894	0.0463862400
	...

Returns:
  - a list reference of genome position / value objects:
    [
      {
        no    					=>  __integer__,          # 1 -> n
        probe_id  			=>	__string__,
        reference_name	=>	__string__,
        position	  		=>	__integer__,
        value						=>	__long__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $probeF  	=   shift;
  my $plot      =   shift;
  my $probes  	=   [];

	if (! -f $probeF) { return $probes }

  open	FILE, "$probeF" or die "No file $probeF $!";
  my $i     		=   0;
  foreach (<FILE>) {
    $i++;
    chomp;
  	my (
  		$probe_id,
  		$reference_name,
  		$position,
  		$value,
  	)		    		=		split (/\s/, $_, 5);

  	$probe_id		=~ 	s/[^\w\-\,]/_/g;
  	$reference_name			=~ s/[^\dxXyY]//;
  	$reference_name			=~ s/^23$/X/;
  	$reference_name			=~ s/^24$/Y/;
  	$position		= 	sprintf "%.0f", $position;	# due to some erroneous .5 in-between pos.
  	$value			=		sprintf "%.4f", $value;

  	if ($reference_name	!~ /^\w\d?$/) 						{ next }
  	if ($position				!~ /^\d{1,9}$/) 					{ next }
  	if ($value					!~ /^\-?\d+?(\.\d+?)?$/) 	{ next }

  	push(
  	  @$probes,
  	  {
  	    no      				=>  $i,
  	    probe_id				=>	$probe_id,
        reference_name  =>	$reference_name,
        position	  		=>	$position,
        value	    			=>	$value,
      }
    );
  }

  close FILE;

  return $probes;

}

################################################################################

sub read_segmentfile {

=pod

Expects:
  - a standard Progenetix segments  file

	sample	chro	start	stop	mean	probes
	GSM481286	1	742429	7883881	-0.1594	699
	GSM481286	1	115673158	115705254	-0.3829	8
	GSM481286	1	115722621	119771659	0.167	424
	GSM481286	1	119776776	162617092	0.4168	1587
	GSM481286	1	162621657	165278686	0.6508	350
	GSM481286	1	165280711	167221337	0.4056	241
	GSM481286	1	167248788	168289603	0.6784	130
	...

Returns:
  - a list reference of genome CNV objects:
    [
      {
        no    					=>  __integer__,    # 1 -> n
        callset_id  		=>	__string__,
        reference_name	=>	__string__,
        start	  				=>	__integer__,
        end	  					=>	__integer__,
        variant_type		=>	__string__,			# DUP, DEL
        info						=>	{
        	value						=>	__long__,
        	svlen						=>	__integer__,
        	probes					=>	__integer__,
        	assembly_id			=>	__string__,			# GRCh36 ...
        	experiment_type	=>	__string__,			# aCGH ...
        },
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $segmentsF =   shift;
  my $plot      =   shift;
  my $segments 	=   [];

	if (! -f $segmentsF) { return $segments }

  open	FILE, "$segmentsF" or die "No file $segmentsF $!";
  my $i     		=   0;
  foreach (<FILE>) {
    $i++;
    chomp;
  	my (
  		$callset_id,
  		$reference_name,
  		$start,
  		$end,
  		$value,
  		$probes,
  	)		    		=		split (/\s/, $_, 6);

  	$callset_id	=~ 	s/[^\w\-\,]/_/g;
  	$reference_name			=~ s/[^\dxXyY]//;
  	$reference_name			=~ s/^23$/X/;
  	$reference_name			=~ s/^24$/Y/;
  	$start			=		sprintf "%.0f", $start;
  	$end  			=		sprintf "%.0f", $end;
  	$probes			=~ 	s/[^\d]//g;
  	$value			=		sprintf "%.4f", $value;

  	if ($reference_name	!~ /^\w\d?$/) 						{ next }
  	if ($start					!~ /^\d{1,9}$/) 					{ next }
  	if ($end						!~ /^\d{1,9}$/) 					{ next }
  	if ($value					!~ /^\-?\d+?(\.\d+?)?$/) 	{ next }

  	push(
  	  @$segments,
  	  {
  	    no      				=>  $i,
  	    callset_id			=>	$callset_id,
        reference_name  =>	$reference_name,
        start	  				=>	1 * $start,
        end		  				=>	1 * $end,
        info						=>	{
        	value					=>	1 * $value,
        	svlen					=>	1 * ($end - $start),
        	probes				=>	1 * $probes,
        },
      }
    );

	}

	for my $i (0..$#{ $segments }) {

    if ($segments->[$i]->{info}->{value} >= $plot->{parameters}->{cna_gain_threshold}) {
      $segments->[$i]->{variant_type}  =   'DUP' }
    elsif ($segments->[$i]->{info}->{value} <= $plot->{parameters}->{cna_loss_threshold}) {
      $segments->[$i]->{variant_type}  =   'DEL' }

  }

	return $segments;

}

################################################################################
########    utility subs    ####    ####    ####    ####    ####    ####    ####
################################################################################

1;
