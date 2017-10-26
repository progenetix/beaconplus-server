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
use PGX::GenomeIntervals::CytobandReader;
use PGX::GenomeIntervals::GenomeIntervals;
use PGX::GenomeIntervals::IntervalStatistics;
use PGX::GenomePlots::HistoPlotter;
use PGX::GenomePlots::PlotParameters;
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
our $submit     =   param('submit');
our $cgi        =   new CGI;

print $cgi->header;
print $cgi->start_html(-title => 'Beacon+ Handover Prototype');

_add_form();
_print_histogram();
_export_callsets();

print $cgi->end_html;





################################################################################
# subs #########################################################################
################################################################################

sub _add_form {

  my $table     =   '
<table>
  <tr><td>User / email</td><td>';
  $table        .=  $cgi->textfield(
                      -name       => 'user',
                      -value      => '',
                      -size       => 30,
                      -maxlength  => 30,
                    );
  $table        .=   '</td></tr>
  <tr><td>Password</td><td>';
  $table        .=  $cgi->password_field(
                      -name       => 'pass',
                      -value      => '',
                      -size       => 30,
                      -maxlength  => 30,
                    );
  $table        .=   '</td></tr>
  <tr><td></td><td>';
    $table      .=  $cgi->submit(
                      -name       => 'submit',
                      -value      => 'Generate Histogram Plot',
                    );
  $table        .=   '</td></tr>
  <tr><td></td><td>';
    $table      .=  $cgi->submit(
                      -name       => 'submit',
                      -value      => 'Export Callset Data',
                    );
  $table        .=  '</td></tr>
</table>';

  print $cgi->start_form;
  print $cgi->hidden(
    -name       => 'accessid',
    -default    => $access_id,
  );
  print $table;
  print $cgi->end_form;

}

################################################################################

sub _print_histogram {

  if ($submit !~ /histo/i) { return }

  my $args        =   {};
  $args->{'-genome'}    =   'grch36';
  $args->{'-binning'}   =   1000000;
  $args->{'-plotid'}    =   'histoplot';
  $args->{'-chr2plot'}  =   '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X';

  $MongoDB::Cursor::timeout = 120000;

  my $tmpcoll   =   MongoDB::MongoClient->new()->get_database( $tempdb )->get_collection($tmpcoll);
  my $tmpdata   =   $tmpcoll->find_one( { _id	=>  $access_id } );

  if ($tmpdata->{query_coll} =~ /_(((?:grch)|(?:hg))\d\d)$/i) {
    $args->{'-genome'}  =   $1 }

  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $tmpdata->{query_db} )->get_collection($tmpdata->{query_coll});
  my $dataQuery =   { $tmpdata->{query_key} => { '$in' => $tmpdata->{query_values} } };
  my $cursor	  =		$datacoll->find( $dataQuery )->fields( { 'info' => 1 } );
  my $data	    =		[ $cursor->all ];

  my $plot      =   new PGX::GenomePlots::Genomeplot($args);
  my $fMaps     =   interval_cnv_frequencies([map{$_->{info}->{statusmaps}} @$data ], $plot->{genomeintervals});
  $plot         =   return_histoplot_svg($plot, $fMaps);

#  print 'Content-type: image/svg+xml'."\n\n";
  print $plot->{svg};

}

################################################################################

sub _export_callsets {

  if ($submit !~ /export/i) { return }

  $MongoDB::Cursor::timeout = 120000;

  my $tmpcoll   =   MongoDB::MongoClient->new()->get_database( $tempdb )->get_collection($tmpcoll);
  my $tmpdata   =   $tmpcoll->find_one( { _id	=>  $access_id } );
  my $datacoll  =   MongoDB::MongoClient->new()->get_database( $tmpdata->{query_db} )->get_collection($tmpdata->{query_coll});
  my $dataQuery =   { $tmpdata->{query_key} => { '$in' => $tmpdata->{query_values} } };
  my $cursor	  =		$datacoll->find( $dataQuery )->fields( { info => 0, _id => 0, updated => 0 } );

  print	JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode([$cursor->all]);

}


1;
