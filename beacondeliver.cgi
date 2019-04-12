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
use BeaconPlus::ConfigLoader;
use BeaconPlus::QueryParameters;
use BeaconPlus::QueryExecution;

use lib './PGX';
use PGX;

=pod

Script for linking callset ids from Beacon+ query to the original data

Please see the associated beaconresponse.md

=cut

use BeaconPlus::ConfigLoader;
use BeaconPlus::QueryParameters;
use BeaconPlus::QueryExecution;

my $config      =   BeaconPlus::ConfigLoader->new();
my $query       =   BeaconPlus::QueryParameters->new($config);

my $error				=		q{};

if ($query->{param}->{accessid}->[0]  =~  /[^\w\-]/) {
  print 'Content-type: text'."\n\n";
  $error				= '¡Wrong access_id parameter '.$query->{param}->{accessid}->[0].'!' }
elsif ($query->{param}->{accessid}->[0]  !~  /\w/) {
  $error				= '¡Missing access_id parameter!' }

if ($error  =~  /\w/) {
  print 'Content-type: text'."\n\n".$error;
  exit;
}

my $pretty     	=   $query->{param}->{jsonpretty}->[0];
if ($pretty !~ /^1|y/) {
  $pretty       =   0 }
else {
  $pretty       =   1 }

################################################################################

our $handover   =    MongoDB::MongoClient->new()->get_database( $config->{handover_db} )->get_collection( $config->{handover_coll} )->find_one( { _id	=>  $query->{param}->{accessid}->[0] } );

# print 'Content-type: text'."\n\n".Dumper($handover->{target_count});
#exit;

_print_histogram();
_export_callsets();
_export_biosamples_individuals();
_export_variants();

exit;

################################################################################
# subs #########################################################################
################################################################################

sub _print_histogram {

  if ($query->{param}->{do}->[0] !~ /cnvhistogram/i) { return }

  $query->{param}->{-chr2plot}  ||=   [1..22];
  $query->{param}->{-plotid}		=		'histoplot';
  $query->{param}->{-text_bottom_left} =  $handover->{source_db}.': '.$handover->{target_count}.' samples';
  
  my $pgx       =   new PGX($query->{param});
  $pgx->pgx_open_handover($config, $query);
  $pgx->pgx_samples_from_accessid($query);
  $pgx->pgx_add_frequencymaps( [ { statusmapsets =>  $pgx->{samples} } ] );
  $pgx->return_histoplot_svg();

  print 'Content-type: image/svg+xml'."\n\n";
  print $pgx->{svg};
  exit;

}

################################################################################

sub _export_callsets {

  if ($query->{param}->{do}->[0] !~ /callsetsdata/i) { return }

  print 'Content-type: application/json'."\n\n";

  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} )->get_collection( $handover->{target_collection} );
  my $dataQuery =   { $handover->{target_key} => { '$in' => $handover->{target_values} } };
  my $cursor	  =		$datacoll->find( $dataQuery )->fields( { _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);
  exit;

}

################################################################################

sub _export_biosamples_individuals {


  if (! grep{ /^$query->{param}->{do}->[0]$/ } qw(biosamplesdata individualsdata)) { return }

  print 'Content-type: application/json'."\n\n";

  my $dataconn  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} );

  my $datacoll  =   $dataconn->get_collection($handover->{target_collection});
  my $cursor	  =		$datacoll->find( { $handover->{target_key} => { '$in' => $handover->{target_values} } } )->fields( { attributes => 0, _id => 0, updated => 0, created => 0 } );

  print	JSON::XS->new->pretty( $pretty )->allow_blessed->convert_blessed->encode([$cursor->all]);
  exit;

}

################################################################################

sub _export_variants {

  if ($query->{param}->{do}->[0] !~ /variants/i) { return }
  print 'Content-type: application/json'."\n\n";

  my $dataconn  =   MongoDB::MongoClient->new()->get_database( $handover->{source_db} );
  my $datacoll  =   $dataconn->get_collection( 'variants' );
  my $cursor    =   {};

  my $key       =   $handover->{target_key};
  my $values    =   $handover->{target_values};
  if ($query->{param}->{do}->[0] =~ /callset/i) { 
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
