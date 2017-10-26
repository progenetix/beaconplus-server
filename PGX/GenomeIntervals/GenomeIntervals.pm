package PGX::GenomeIntervals::GenomeIntervals;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(make_genome_intervals get_reference_base_limits get_genome_basecount);

################################################################################

sub make_genome_intervals {

=pod

Expects:
  - a list reference of genome interval objects, usually representing cytobands:
    [
      {
        no    =>  __integer__,            # not used
        reference_name  =>	__string__,
        start	=>	__integer__,
        end	  =>	__integer__,
        stain	=>	__string__,            # not used
        label	=>	__string__,            # not used
      },
      {
      ...
      },
    ]
  - the binning size in bases (optional; defaults to 1000000)

Returns:
  - a list of genome intervals in a similar structure, based on the binning size
  - the intervals are per reference, i.e. each reference starts with a new interval,
    leading to the last interval < binning size
    [
      {
        no    =>  __integer__,          # 1 -> n
        reference_name  =>	__string__,
        start	=>	__integer__,
        end	  =>	__integer__,
        label	=>	__string__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $intSize;
  ($plotPars->{CYTOBANDS}, $intSize)     =   @_;

  if ($intSize !~ /^\d{3,9}$/) { $intSize = 1000000 }

  my $refLims   =   get_reference_base_limits($plotPars->{CYTOBANDS});
  my $gi        =   [];

  # references are sorted with numerical ones first, then others (e.g. 1 -> 22, X, Y)
  my @refNames  =   ((sort {$a <=> $b } grep{ /^\d\d?$/ } keys %$refLims), (sort grep{ ! /\d/ } keys %$refLims));

  my $intI      =   1;
  for my $i (0..$#refNames) {

    my $refName =   $refNames[$i];
    my $start   =   $refLims->{ $refName }->[0];
    my $end     =   $intSize - 1;

    while ($start <= $refLims->{ $refName }->[1]) {

      # adjusting the end of the last interval
      if ($end > $refLims->{ $refName }->[1]) { $end = $refLims->{ $refName }->[1] };

      push(
        @$gi,
        {
          no    =>  $intI,
          reference_name  =>  $refName,
          start =>  $start,
          end   =>  $end,
          label =>  $refName.':'.$start.'-'.$end,
        }
      );

      $start    +=  $intSize;
      $end      +=  $intSize;
      $intI++;

  }}

  return $gi;

}

################################################################################

sub get_reference_base_limits {

=pod

Expects:
  - a list reference of genome interval objects:
    [
      {
        reference_name  =>  __string__,
        start           =>  __integer__,
        end             =>  __integer__,
        (others)...
      },
      {
        ...
      },
    ]

Returns:
  - a hash reference with each value consisting of a list of 2 integers, representing
    the reference's min and max bases:
    [
      { __string__  =>  [  __integer__,  __integer__ ] },
      { ... },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $refLims   =   {};
  my $allRefs   =   $_[0];

  foreach my $ref (map{ $_->{reference_name} } @$allRefs) {
    my @refRefs =   grep{ $_->{reference_name} =~ /^$ref$/i } @$allRefs;
    my @bases   =   sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
    $refLims->{$ref}  =   [ $bases[0], $bases[-1] ];
  }

  return $refLims;

}

################################################################################

sub get_genome_basecount {

  my $genBases;
  my ($allRefs, $chr2plot)   =   @_;
#   my %refNames  =   map { $_->{reference_name} => 1 } @$allRefs;
  my %refNames  =   map { $_ => 1 } @$chr2plot;

  foreach my $ref (keys %refNames) {
    my @refRefs =   grep{ $_->{reference_name} =~ /^$ref$/i } @$allRefs;
    my @bases   =   sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
    $genBases   +=   ($bases[-1] - $bases[0]);
  }

  return $genBases;

}


################################################################################
########    utility subs    ####    ####    ####    ####    ####    ####    ####
################################################################################

sub _genome_names_to_grch {
  my %genome_names   =   (
    hg18    =>  'grch36',
    hg19    =>  'grch37',
    hg38    =>  'grch38',
  );
  if ($_[0] =~ /^grch\d\d$/i) {
    return lc($_[0]) }
  else {
    return $genome_names{lc($_[0])} }
}

################################################################################

sub _genome_names_to_hg {
  my %genome_names   =   (
    grch36  =>  'hg18',
    grch37  =>  'hg19',
    grch38  =>  'hg38',
  );
  if ($_[0] =~ /^hg\d\d$/i) {
    return lc($_[0]) }
  else {
    return $genome_names{lc($_[0])} }
}

1;
