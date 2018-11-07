#!/usr/bin/perl

# © 2000-2018 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(param);

use File::Basename;
use JSON;
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;

use YAML::XS qw(LoadFile);

my $here_path   =   File::Basename::dirname( eval { ( caller() )[1] } );
our $config     =   LoadFile($here_path.'/rsrc/config.yaml') or die print 'Content-type: text'."\n\n¡No config.yaml file in this path!";

my $querytype   =   param('querytype');
my $datasets    =   [ param('dataset_ids') ];

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

my %ontologyIds;
my $biosQ       =   {};
my $querytext;
my $queryregex;

if (param( 'querytext')) {
  $querytext    =   param( 'querytext');
  $queryregex   =   qr/$querytext/i;
  $biosQ        =   { "id" => { '$regex' => $queryregex } } }

#print Dumper($biosQ);

foreach my $datasetId ( @{ $config->{ dataset_names }}) {

  if (! grep{ /$datasetId/ } @dbList) { next }

  my $dbconn    =   $dbClient->get_database( $datasetId );

  my $datasetI  =       {
    identifier  => 'org.progenetix:'.$config->{ beacon_id }.':'.$datasetId,
    datasetId   => $datasetId,
    name        => $datasetId,
  };

  foreach my $coll (@{ $config->{ collection_names } }) {
    $datasetI->{ $coll.'_count' } =   $dbconn->get_collection($coll)->count;
  }

  ##############################################################################

  my $cursor    =   $dbconn->get_collection('biosubsets')->find( $biosQ )->fields({ id => 1, label => 1, count => 1});
  foreach ($cursor->all) {
    $ontologyIds{ $_->{id} } = $_->{id}.': '.$_->{label};
    if (@{ $config->{ dataset_names }} == 1) {
      $ontologyIds{ $_->{id} } .= ' ('.$_->{count}.')';
    }
  }

  $datasetI->{ info }->{ ontology_terms } =   [ map{ { $_ => $ontologyIds{ $_ } } } sort keys %ontologyIds ];

  ##############################################################################

  push(
    @{ $beaconInfo->{dataset} },
    $datasetI
  );

}

$beaconInfo->{supportedRefs}  =   $config->{ genome_assemblies };

if ($querytype =~ /ontologyid/) {
  $beaconInfo  =   [ map{ { id => $_, infolabel => $ontologyIds{$_} } } sort keys %ontologyIds ] }

print JSON::XS->new->pretty( 1 )->allow_blessed->convert_blessed->encode($beaconInfo)."\n";

exit;

1;
