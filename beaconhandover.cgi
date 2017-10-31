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
# use PGX::GenomeIntervals::CytobandReader;
# use PGX::GenomeIntervals::GenomeIntervals;
# use PGX::GenomeIntervals::IntervalStatistics;
# use PGX::GenomePlots::HistoPlotter;
# use PGX::GenomePlots::PlotParameters;
# use PGX::GenomePlots::Genomeplot;
#
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

  print $cgi->start_form(-action => 'beacondeliver.cgi');
  print $cgi->hidden(
    -name       => 'accessid',
    -default    => $access_id,
  );
  print $table;
  print $cgi->end_form;

}

################################################################################

1;
