package PGX::GenomePlots::Genomeplot;

use Data::Dumper;
use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomePlots::PlotParameters;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(new);

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



1;
