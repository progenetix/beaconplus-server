#!/usr/bin/perl -w

# Progenetix & arrayMap site scripts
# © 2000-2018 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);
use JSON::XS;
use MongoDB;
use MongoDB::MongoClient;

my $request	    =		param('querytext');
my $db					=		'progenetix';
my $collection  =		'genespans';
my $minChars		=		param('minchars');
my $debug				=		param('debug');

if ($minChars !~ /^\d+$/) { $minChars = 1 }

if ($debug > 0) {
	print 'Content-type: text/plain'."\n\n" }

if ($request !~ /\w{$minChars,127}/ || $collection !~ /\w/i) {
	print 'Content-type: text/plain'."\n\n";
	exit;
}

################################################################################

my $cursor			=		MongoDB::MongoClient->new()->get_database( $db )->get_collection( $collection )->find( { "gene_symbol" => qr/^$request/i } )->fields( {_id => 0, gene_symbol => 1, reference_name => 1, cds_start_min => 1, cds_end_max => 1, } );

print	'Content-type: application/json'."\n\n";

print	param('callback').'({"genes":['.join(',', (map{ JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode($_) } $cursor->all ))."]});\n";

exit;

1;
