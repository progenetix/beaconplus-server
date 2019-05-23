#!/usr/bin/perl

# © 2000-2019 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);

use File::Basename;
use JSON;
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;

use BeaconPlus::ConfigLoader;
my $config      =   BeaconPlus::ConfigLoader->new();
my $querytype   =   param('querytype');
my $autoc       =   param('autocomplete');
my $datasets    =   [ param('datasetIds') ];

if ($datasets->[0] =~ /\w\w\w/) {
  $config->{ dataset_names }  =   [];
  foreach (grep{ /\w\w\w/ } @$datasets) {
    push(@{ $config->{ dataset_names } }, $_);
  }
}

if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

if ($querytype =~ /get_datasetids/) {

  print JSON::XS->new->pretty( 1 )->allow_blessed->convert_blessed->encode($config->{ dataset_names })."\n";

  exit;

}

=pod

Please see the associated beaconresponse.md

=cut

my $dbClient    =   MongoDB::MongoClient->new();
my @dbList      =   $dbClient->database_names();

my $beaconInfo  =   {
  beaconId      =>  $config->{ beacon_id },
  provider      =>  $config->{ provider },
  description   =>  $config->{ description },
  url           =>  $config->{ url },
  sameAs        =>  $config->{ url_alt },
  logo          =>  $config->{ url_logo },
  potentialActions   =>  $config->{ actions },
  dataset       =>  [],
};

my $ontologyIds;
my $biosQ       =   {};
my $querytext;
my $queryregex;

if (param( 'querytext')) {
  $querytext    =   param( 'querytext');
  $queryregex   =   qr/$querytext/i;
  $biosQ        =   { '$or' =>  [
                      { "id"    => { '$regex' => $queryregex } },
                      { "label" => { '$regex' => $queryregex } },
                    ] };
} else {
  print JSON::XS->new->pretty( 1 )->allow_blessed->convert_blessed->encode($beaconInfo)."\n";
  exit;
}

#print Dumper($biosQ);

foreach my $datasetId ( @{ $config->{ dataset_names }}) {

  if (! grep{ /$datasetId/ } @dbList) { next }

  my $dbconn    =   $dbClient->get_database( $datasetId );

  my $datasetI  =       {
    identifier  => 'org.progenetix:'.$config->{ beacon_id }.':'.$datasetId,
    datasetId   => $datasetId,
    name        => $datasetId,
  };

  foreach my $coll (values %{ $config->{ collection_names } }) {
    $datasetI->{ $coll.'_count' } =   $dbconn->get_collection($coll)->count;
  }

  ##############################################################################
  
  my $collName  =   'biosubsets';
  if ($querytype  =~ /referenceid/) {
    $collName  =   'datacollections' }

  my $cursor    =   $dbconn->get_collection($collName)->find( $biosQ )->fields({ id => 1, label => 1, child_terms => 1, count => 1, _id => 0});
  my @subsets   =   $cursor->all;
  foreach my $sub (@subsets) {
    $sub->{label_short}   =  $sub->{label};
    $sub->{label_short}   =~ s/^(.{20,45}?)[\s\,].*?$/$1.../;
    if ($sub->{id} =~ /\+$/) {
    	$sub->{child_terms}	=	join(',', grep{ $_ !~ /\+/ } @{ $sub->{child_terms} }) }
   	else {
   		$sub->{child_terms} =  $sub->{id} }
    $ontologyIds->{ $sub->{id} }  =   $sub;
  }

  $datasetI->{ info }->{ ontology_terms } =   [ map{ { $_ => $ontologyIds->{ $_ } } } sort keys %{ $ontologyIds } ];

  ##############################################################################

  push(
    @{ $beaconInfo->{dataset} },
    $datasetI
  );

}

$beaconInfo->{supportedRefs}  =   $config->{ genome_assemblies };

if ($querytype =~ /ontologyid|referenceid/) {
  $beaconInfo  =   [ map{ $ontologyIds->{$_} } (grep{ /[^\-\+]/ } sort keys %{ $ontologyIds } )] }

if ($autoc =~ /1|y/i) {
  print param('callback').'({"data":['.join(',', (map{ JSON::XS->new->pretty( 0 )->allow_blessed->convert_blessed->encode($_) } @$beaconInfo )).']});'."\n" }
else {
  print JSON::XS->new->pretty( 1 )->allow_blessed->convert_blessed->encode($beaconInfo)."\n" }

exit;

1;
