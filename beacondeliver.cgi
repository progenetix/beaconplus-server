#!/usr/bin/perl

# Beacon+ support scripts
# © 2017 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(:standard param *table);

use JSON;
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;
use UUID::Tiny;

# local packages
BEGIN {
 push @INC, '/Library/WebServer/cgi-bin';
}
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::HistoPlotter;
use PGX::GenomePlots::Genomeplot;

=pod

Script for linking callset ids from Beacon+ query to the original data

Please see the associated beaconresponse.md

=cut

our $tempdb     =   'progenetix';
our $tmpcoll    =   'querybuffer';

#if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

# parameters
our $access_id  =   param('accessid');
our $dataStyle  =   param('datastyle');
our $todo       =   param('do');
our $genome     =   param('assembly_id');
our $chr2plot   =   param('chr2plot');
our $cgi        =   new CGI;

_print_histogram();
_export_callsets();
_export_biosamples();

################################################################################
# subs #########################################################################
################################################################################

sub _print_histogram {

  if ($todo !~ /histo/i) { return }


  my $args      =   {};
  $args->{'-genome'}    ||= $genome     ||= 'grch36';
  $args->{'-chr2plot'}  ||= $chr2plot   ||= join(',', 1..22,'X');
  $args->{'-binning'}   =   1000000;
  $args->{'-plotid'}    =   'histoplot';
  $args->{'-plottype'}       =   'histogram';

  $MongoDB::Cursor::timeout = 120000;

  my $tmpcoll   =   MongoDB::MongoClient->new()->get_database( $tempdb )->get_collection($tmpcoll);
  my $tmpdata   =   $tmpcoll->find_one( { _id =>  $access_id } );
  if ($tmpdata->{query_coll} =~ /_(((?:grch)|(?:hg))\d\d)$/i) {
    $args->{'-genome'}  =   $1 }

  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $tmpdata->{query_db} )->get_collection($tmpdata->{query_coll});
  my $dataQuery =   { $tmpdata->{query_key} => { '$in' => $tmpdata->{query_values} } };
  my $cursor    =   $datacoll->find( $dataQuery )->fields( { info => 1, statusmaps => 1 } );
  my $callsets  =   [ $cursor->all ];

  $args->{'-text_bottom_left'}  =   scalar(@$callsets).' samples';

  if ($dataStyle =~ /old/i) {
    $callsets   =   [ map{  { info => { statusmaps => $_->{statusmaps}}} } @$callsets ];
  }

  my $plot      =   new PGX::GenomePlots::Genomeplot($args);
  $plot->plot_add_frequencymaps($callsets);
  $plot->return_histoplot_svg();

  print 'Content-type: image/svg+xml'."\n\n";
  print $plot->{svg};

}

################################################################################

sub _export_callsets {

  if ($todo !~ /callsets/i) { return }

  print 'Content-type: application/json'."\n\n";

  $MongoDB::Cursor::timeout = 120000;

  my $tmpcoll   =   MongoDB::MongoClient->new()->get_database( $tempdb )->get_collection($tmpcoll);
  my $tmpdata   =   $tmpcoll->find_one( { _id =>  $access_id } );
  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $tmpdata->{query_db} )->get_collection($tmpdata->{query_coll});
  my $dataQuery =   { $tmpdata->{query_key} => { '$in' => $tmpdata->{query_values} } };
  my $cursor    =   $datacoll->find( $dataQuery )->fields( { info => 0, _id => 0, updated => 0, created => 0 } );

  print JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode([$cursor->all]);

}

################################################################################

sub _export_biosamples {

  if ($todo !~ /biosamples/i) { return }

  print 'Content-type: application/json'."\n\n";

  $MongoDB::Cursor::timeout = 120000;

  my $tmpcoll   =   MongoDB::MongoClient->new()->get_database( $tempdb )->get_collection($tmpcoll);
  my $tmpdata   =   $tmpcoll->find_one( { _id =>  $access_id } );
  my $dataconn  =   MongoDB::MongoClient->new()->get_database( $tmpdata->{query_db} );
  my $datacall  =   $dataconn->run_command([
                      "distinct"=>  $tmpdata->{query_coll},
                      "key"     =>  'biosample_id',
                      "query"   =>  { $tmpdata->{query_key} => { '$in' => $tmpdata->{query_values} } },
                    ]);
  my $biosids   =   $datacall->{values};
  my $datacoll  =   $dataconn->get_collection('biosamples');
  my $cursor    =   $datacoll->find( { id => { '$in' => $biosids } } )->fields( { attributes => 0, _id => 0, updated => 0, created => 0 } );

  print JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode([$cursor->all]);

}


1;
