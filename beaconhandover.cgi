#!/usr/bin/perl

# Beacon+ support scripts
# © 2017 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);

use JSON;
use List::Util qw(min max);
use List::MoreUtils qw(any apply);
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;
use UUID::Tiny;

=pod

Script for linking callset ids from Beacon+ query to the original data

Please see the associated beaconresponse.md

=cut

if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

my $args        =   {};

my $tempdb      =   'progenetix';
my $tmpcoll     =   'querybuffer';
my $access_id   =   param('accessid');

$MongoDB::Cursor::timeout = 120000;

my  $tmpcoll    =   MongoDB::MongoClient->new()->get_database( $tempdb )->get_collection($tmpcoll);

my $tmpdata     =   $tmpcoll->find_one( { _id	=>  $access_id } );

my $datacoll    =   MongoDB::MongoClient->new()->get_database( $tmpdata->{query_db} )->get_collection($tmpdata->{query_coll});

my $dataQuery   =   { $tmpdata->{query_key} => { '$in' => $tmpdata->{query_values} } };

my $cursor	    =		$datacoll->find( $dataQuery )->fields( { 'statusmaps.dupmap' => 1 } );
my $data	      =		[ $cursor->all ];
print Dumper($data);







1;
