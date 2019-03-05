#!/usr/bin/perl

# Beacon+ support scripts
# © 2017-20018 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);

use JSON;
use MongoDB;
use MongoDB::MongoClient;
$MongoDB::Cursor::timeout = 120000;

use Data::Dumper;

# local packages
use beaconPlus::ConfigLoader;
use beaconPlus::QueryParameters;
use beaconPlus::QueryExecution;
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::HistoPlotter;
use PGX::GenomePlots::Genomeplot;

=pod

Script for linking callset ids from Beacon+ query to the original data

Please see the associated beaconresponse.md

=cut

use beaconPlus::ConfigLoader;
use beaconPlus::QueryParameters;
use beaconPlus::QueryExecution;

my $config      =   beaconPlus::ConfigLoader->new();
my $query       =   beaconPlus::QueryParameters->new($config);

# parameters
my $access_id   =   param('accessid');
our $todo       =   param('do');
our $pretty     =   param('jsonpretty');

if ($access_id  =~  /[^\w\-]/) {
  print 'Content-type: text'."\n\n";
  print '¡Wrong access_id parameter '.$access_id.'!';
  exit;
}

if ($access_id  !~  /\w/) {
  print 'Content-type: text'."\n\n";
  print '¡Missing access_id parameter!';
  exit;
}

if ($pretty !~ /^1|y/) {
  $pretty       =   0 }

our $handover   =    MongoDB::MongoClient->new()->get_database( $config->{handover_db} )->get_collection( $config->{handover_coll} )->find_one( { _id	=>  $access_id } );
#print 'Content-type: text'."\n\n";
# print Dumper($todo);
# print Dumper($handover);
#exit;
_print_histogram();
_export_callsets();
_export_biosamples_individuals();
_export_variants();
print 'Content-type: text'."\n\n";
exit;


################################################################################
# subs #########################################################################
################################################################################

sub _print_histogram {

  if ($todo !~ /cnvhistogram/i) { return }

  my $plotargs    =   {
    -chr2plot     =>   join(',', 1..22),
    -plotid       =>  'histoplot',
    -do_plottype  =>  'histogram',
    -size_plotimage_w_px  =>  800,
    -size_plotarea_h_px   =>  100,
    -size_chromosome_w_px =>  12,
    -size_text_px =>  9,
    -size_title_left_px   =>   180,
    -size_text_title_left_px  =>  10,
  };

  foreach my $par (grep{ /^\-\w/ } keys %{ $query->{param} }) {
    if (grep{$par eq $_ } qw(-chr2plot -markers)) {
      $plotargs->{$par} =   join(',', @{ $query->{param}->{$par} }) }
    else {
      $plotargs->{$par} =   $query->{param}->{$par}->[0] }
  }


  my $pgx         =   new PGX::GenomePlots::Genomeplot($plotargs);
  $pgx->pgx_open_handover($config, $query);
  $pgx->pgx_samples_from_accessid($query);

  my $thisargs    =   { map{ $_ => $plotargs->{$_} } keys %{$plotargs} };
  delete $thisargs->{-size_title_left_px};
  $thisargs->{-text_bottom_left}    =  $pgx->{handover}->{source_db}.': '.scalar(@{ $pgx->{samples} }).' samples';

  my $plot        =   new PGX::GenomePlots::Genomeplot($thisargs);
  $plot->pgx_add_frequencymaps( [ { statusmapsets =>  $pgx->{samples} } ] );
  $plot->return_histoplot_svg();


  print 'Content-type: image/svg+xml'."\n\n";
  print $plot->{svg};
  exit;

}

################################################################################

sub _export_callsets {

  if ($todo !~ /callsetsdata/i) { return }

  print 'Content-type: application/json'."\n\n";

  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} )->get_collection( $handover->{target_collection} );
  my $dataQuery =   { $handover->{target_key} => { '$in' => $handover->{target_values} } };
  my $cursor	  =		$datacoll->find( $dataQuery )->fields( { _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);
  exit;

}

################################################################################

sub _export_biosamples_individuals {


  if (! grep{ /^$todo$/ } qw(biosamplesdata individualsdata)) { return }

  print 'Content-type: application/json'."\n\n";

  my $dataconn  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} );

  my $datacoll  =   $dataconn->get_collection($handover->{target_collection});
  my $cursor	  =		$datacoll->find( { $handover->{target_key} => { '$in' => $handover->{target_values} } } )->fields( { attributes => 0, _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);
  exit;

}

################################################################################

sub _export_variants {

  if ($todo !~ /variants/i) { return }
  print 'Content-type: application/json'."\n\n";

  my $dataconn  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} );
  my $datacoll  =   $dataconn->get_collection( 'variants' );
  my $cursor    =   {};

  my $key       =   $handover->{target_key};
  my $values    =   $handover->{target_values};
  if ($todo =~ /callset/i) { 
    $key        =   'callset_id';
    my $distincts =   $dataconn->run_command([
                      "distinct"=>  'callsets',
                      "key"     =>  'id',
                      "query"   =>  { _id => { '$in' => $handover->{target_values} } },
                    ]);

    $values     =   $distincts->{values};

  }
  $cursor	      =		$datacoll->find( { $key => { '$in' => $values } } )->fields( { _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);
  exit;

}


1;
