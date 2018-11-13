#!/usr/bin/perl

# Beacon+ support scripts
# © 2017-20018 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(:standard param *table);

use JSON;
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;
use YAML::XS qw(LoadFile);

# local packages
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::HistoPlotter;
use PGX::GenomePlots::Genomeplot;

=pod

Script for linking callset ids from Beacon+ query to the original data

Please see the associated beaconresponse.md

=cut

my $here_path   =   File::Basename::dirname( eval { ( caller() )[1] } );
our $config     =   LoadFile($here_path.'/rsrc/config.yaml') or die print 'Content-type: text'."\n\n¡No config.yaml file in this path!";

#print 'Content-type: text'."\n\n";

# parameters
my $access_id   =   param('accessid');
our $todo       =   param('do');
our $pretty     =   param('jsonpretty');
our $cgi        =   new CGI;

if ($access_id  =~  /[^\w\-]/) {

  print 'Content-type: text'."\n\n";
  print 'Wrong of missing access_id parameter '.$access_id;

  exit;
}

if ($pretty !~ /^1|y/) {
  $pretty       =   0 }

$MongoDB::Cursor::timeout = 120000;

our $handover   =    MongoDB::MongoClient->new()->get_database( $config->{handover_db} )->get_collection( $config->{handover_coll} )->find_one( { _id	=>  $access_id } );
#print Dumper($handover);
#exit;
_print_histogram();
_export_callsets();
_export_biosamples();
_export_variants();

################################################################################
# subs #########################################################################
################################################################################

sub _print_histogram {

  if ($todo !~ /histo/i) { return }

  my $args      =   {};
  $args->{'-plotid'}    =   'histoplot';
  $args->{'-do_plottype'}   =   'histogram';

  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} )->get_collection( $handover->{target_collection} );
  my $dataQuery =   { $handover->{target_key} => { '$in' => $handover->{target_values} } };
  my $cursor	  =		$datacoll->find( $dataQuery )->fields( { 'info' => 1 } );
  my $callsets	=		[ $cursor->all ];

  $callsets     =   [ grep{ exists $_->{info}->{statusmaps} } @$callsets ];
  $args->{'-text_bottom_left'}  =   scalar(@$callsets).' samples';

  # print 'Content-type: text'."\n\n";
  # print Dumper($callsets);

  my $plot      =   new PGX::GenomePlots::Genomeplot($args);
  $plot->plot_add_frequencymaps( [ { statusmapsets =>  $callsets } ] );
  $plot->return_histoplot_svg();

  print 'Content-type: image/svg+xml'."\n\n";
  print $plot->{svg};

}

################################################################################

sub _export_callsets {

  if ($todo !~ /callsets/i) { return }

  print 'Content-type: application/json'."\n\n";

  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} )->get_collection( $handover->{target_collection} );
  my $dataQuery =   { 'id' => { '$in' => $handover->{query_values} } };
  my $cursor	  =		$datacoll->find( $dataQuery )->fields( { _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);

}

################################################################################

sub _export_biosamples {

  if ($todo !~ /biosamples/i) { return }

  print 'Content-type: application/json'."\n\n";

  my $dataconn  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} );
  
  my $datacoll  =   $dataconn->get_collection($handover->{target_collection});
  my $cursor	  =		$datacoll->find( { $handover->{target_key} => { '$in' => $handover->{target_values} } } )->fields( { attributes => 0, _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);

}

################################################################################

sub _export_variants {

  if ($todo !~ /variants/i) { return }
  print 'Content-type: application/json'."\n\n";

  my $dataconn  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} );
  my $datacoll  =   $dataconn->get_collection( 'variants' );
  my $cursor    =   {};

  $cursor	      =		$datacoll->find( { $handover->{target_key} => { '$in' => $handover->{target_values} } } )->fields( { _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);

}


1;
