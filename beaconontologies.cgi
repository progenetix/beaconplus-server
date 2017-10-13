#!/usr/bin/perl -w

# Progenetix & arrayMap site scripts
# © 2000-2017 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);
#use Data::Dumper qw(Dumper);
use JSON;
use MongoDB;
use MongoDB::MongoClient;

my $datasetId		=		param('dataset_id');
my $collection  =		'bioontologies';
my $afqfield	  =		param('afqfield');
my $request	    =		param('querytext');

if ($datasetId !~ /_ga4gh$/ ) 	{ $datasetId .= '_ga4gh' }
if ($afqfield !~ /\w\w\w/ ) { $afqfield = 'term_id' }
if ($request !~ /^[\w\_\:]{1,64}/ ) { $request = '..' }

################################################################################

my $dbconn			=		MongoDB::MongoClient->new()->get_database( $datasetId );
my $fields		  =		{ _id	=>	0, term_id => 1, infolabel => 1 };
my $cursor		  =		$dbconn->get_collection( $collection )->find( { $afqfield => qr/$request/i } )->fields( $fields );
my @items       =   $cursor->all;

@items  =   sort { $a->{term_id} cmp $b->{term_id} } @items;

print	'Content-type: application/json'."\n\n";
print '['.join(',', (map{ JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode($_) } @items ))."];\n";

exit;

1;
