package PGX::GenomePlots::Genomeplot;

use Data::Dumper;
use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomePlots::PlotParameters;
use PGX::FileUtilities::ArrayfileReader;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(new plot_add_probedata plot_add_segmentdata plot_add_probedata_fracb plot_add_segmentdata_fracb);

################################################################################

sub new {

  my $class     =   shift;
  my $args      =   shift;
  my $self      =   {
    parameters  =>  args_modify_plot_parameters(read_plot_defaults(), $args),
    cytobands   =>  read_cytobands($args->{'-genome'}),
    plotid      =>  ($args->{'-plotid'} !~ /^\w+?/ ? 'genomeplot' : $args->{'-plotid'}),
    svg         =>  q{},
    Y           =>  0,
  };
  bless $self, $class;
  $self->{genomeintervals}      =   make_genome_intervals(
                                      $self->{cytobands},
                                      $args->{'-binning'},
                                    );
  $self->{referencebounds}      =   get_reference_base_limits($self->{cytobands});
  $self->{genomesize}           =   get_genome_basecount(
                                      $self->{cytobands},
                                      $self->{parameters}->{chr2plot},
                                    );

  return $self;

}

################################################################################

sub plot_add_probedata {

  my $probefile =   shift;
  my $plot      =   shift;
  $plot->{probedata}    =   read_probefile($probefile, $plot);
  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_add_segmentdata {

  my $segfile   =   shift;
  my $plot      =   shift;
  $plot->{segmentdata}  =   read_segmentfile($segfile, $plot);
  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_add_probedata_fracb {

  my $probefile =   shift;
  my $plot      =   shift;
  $plot->{probedata_fracb}    =   read_probefile($probefile, $plot);
  return $plot;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub plot_add_segmentdata_fracb {

  my $segfile   =   shift;
  my $plot      =   shift;
  $plot->{segmentdata_fracb}  =   read_segmentfile($segfile, $plot);
  return $plot;

}

1;
