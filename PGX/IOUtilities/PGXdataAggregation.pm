package PGX::IOUtilities::PGXdataAggregation;

use Data::Dumper;
use Math::Random qw(random_normal);
use PGX::GenomePlots::PlotParameters;
use PGX::Helpers::UtilityLibs;

require Exporter;
@ISA    =   qw(Exporter);
@EXPORT =   qw(
  pgx_open_handover
  pgx_samples_from_accessid
  pgx_add_segments_from_variants
  pgx_create_samples_from_segments
  pgx_callset_labels_from_biosamples
  pgx_callset_labels_from_file
  pgx_create_sample_collections
);

################################################################################

sub pgx_open_handover {

  my $pgx       =   shift;
  my $config    =   shift;
  my $query     =   shift;
  
  if ($query->{param}->{accessid}->[0] !~ /..../) { return $pgx }

  $pgx->{handover}  =    MongoDB::MongoClient->new()->get_database( $config->{handover_db} )->get_collection( $config->{handover_coll} )->find_one( { _id	=>  $query->{param}->{accessid}->[0] } );
  $pgx->{dataconn}  =   MongoDB::MongoClient->new()->get_database( $pgx->{handover}->{source_db} );

 return $pgx;

}

################################################################################

sub pgx_samples_from_accessid {

  my $pgx       =   shift;
  my $query     =   shift;
  if (! $pgx->{handover}) { return $pgx }
  if ($pgx->{handover}->{target_collection} ne 'callsets') { return $pgx }

  my $cscoll      =   $pgx->{dataconn}->get_collection( $pgx->{handover}->{target_collection} );
  my $dataQuery   =   { $pgx->{handover}->{target_key} => { '$in' => $pgx->{handover}->{target_values} } };
  my $cursor	    =		$cscoll->find( $dataQuery )->fields( { _id => 1, id => 1, biosample_id => 1, 'info.statusmaps.dupmap' => 1, 'info.statusmaps.delmap' => 1 } );
  my $callsets	  =		[ $cursor->all ];
  $callsets       =   [ grep{ exists $_->{info}->{statusmaps} } @$callsets ];
  if ($query->{param}->{'-randno'}->[0] > 0) {
    $callsets     =   RandArr($callsets, $query->{param}->{'-randno'}->[0]) }
 
  $pgx->{samples} =   [ 
    map{
      { 
        id            => $_->{id}, 
        biosample_id  => $_->{biosample_id}, 
        statusmaps    => $_->{info}->{statusmaps}, 
      }
    } @$callsets
  ];
 
  return $pgx;

}

################################################################################

sub pgx_create_samples_from_segments {

  my $pgx       =   shift;
  
  if (! $pgx->{segmentdata}) { return $pgx }
  
  my %csIds     =   map{ $_->{callset_id} => 1 } @{ $pgx->{segmentdata} };

  $pgx->{samples}   ||= [];
 
  foreach my $csId (keys %csIds) {

    # a bit awkward re-assignment of segmentsdata
    my $segments    =   [ grep{ $_->{callset_id} eq $csId } @{ $pgx->{segmentdata} } ];
    $pgx->plot_segments_add_statusmap($segments);
    push(
      @{ $pgx->{samples} },
      {
        id        =>  $csId,
        statusmaps    =>  $pgx->{statusmaps},
        variants  =>  $segments,
        name      =>  $csId,
      }
    );

  }

  return $pgx;
  
}

################################################################################

sub pgx_add_segments_from_variants {

  my $pgx       =   shift;

  if (! $pgx->{dataconn}) { return $pgx }

  my %csIds     =   map{ $_->{id} => 1 } @{ $pgx->{samples} };
  my $cursor    =   $pgx->{dataconn}->get_collection( 'variants' )->find( { callset_id => { '$in' => [ keys %csIds ] } } )->fields( { _id => 0, updated => 0, created => 0, digest => 0, reference_bases => 0, mate_name => 0} );
  my $vars      =   [ $cursor->all ];

  for my $i (0..$#{ $pgx->{samples} }) {
 
    $pgx->{samples}->[$i]->{variants} =   [ grep{ $_->{callset_id} eq $pgx->{samples}->[$i]->{id} } @$vars ];

  }
  
  return $pgx;

}

################################################################################

sub pgx_create_sample_collections {

  my $pgx       =   shift;
  
  $pgx->{samplecollections} =   [];  
  
  my %sortKeys  =   map{ $_->{sortkey} => 1 } @{ $pgx->{samples} };

  # creation of the groups

  if (scalar(keys %sortKeys) == 1) {

    my $groupTag  =   scalar(@{ $pgx->{samples} }).' samples';
    if ((keys %sortKeys)[0] =~ /.../) {
      $groupTag   =   (keys %sortKeys)[0].' ('.$groupTag.')';
    }
    $pgx->{samplecollections} =   [ {
      labels        =>  [ q{} ],
      name          =>  $groupTag,
      statusmapsets =>  [ map{  { statusmaps => $_->{statusmaps} } } @{ $pgx->{samples} } ],
    } ];

  } else {

    foreach my $sortKey (keys %sortKeys) {

      my @theirIndex  =   grep{ $pgx->{samples}->[$_]->{sortkey} eq $sortKey } 0..$#{ $pgx->{samples} };

      my $labelColor  =   random_hexcolor();
      if ( @theirIndex < $pgx->{parameters}->{min_group_no} ) {
        $labelColor  =   '#cccccc' }
      my $label     =   {
        label_text  =>  $sortKey,
        label_link  =>  q{},
        label_color =>  $labelColor,
      };
      $sortKey  =~  s/^.+?\:\:?//;
      foreach (@theirIndex) {
        $pgx->{samples}->[$_]->{labels} =  [ $label ];
        $pgx->{samples}->[$_]->{name}   =  $pgx->{samples}->[$_]->{id}.' ('.$sortKey.')';
        $pgx->{samples}->[$_]->{name}   =~  s/^.+?\:\://g;
        $pgx->{samples}->[$_]->{name}   =~  s/^.+?\:\://g;
      }

      if ( @theirIndex < $pgx->{parameters}->{min_group_no} ) { next }

      my $theseSamples =  [@{ $pgx->{samples} }[@theirIndex]];
      my $name      =   $theseSamples->[0]->{sortlabel}.' ('.scalar(@$theseSamples).')';

      push(
        @{$pgx->{samplecollections}},
        {
          labels        =>  [ $label ],
          name          =>  $name,
          statusmapsets =>  [ map{  { statusmaps => $_->{statusmaps} } } @{ $theseSamples } ],
        },
      );

    }
  }
  
  return $pgx;

}

################################################################################

sub pgx_callset_labels_from_biosamples {

  my $pgx       =   shift;
  my $query     =   shift;
  
  if (! $pgx->{dataconn}) { return $pgx }

  my ($groupAttr, $groupType) =   split('::', $query->{param}->{-grouping}->[0]);
  if ($groupAttr !~ /.../)  { $groupAttr = 'biocharacteristics' }
  if ($groupType !~ /.../)  { $groupType = 'xxxx' }
  my %biosIds   =   map{ $_->{biosample_id} => 1 } @{ $pgx->{samples} };
  my $bscoll    =   $pgx->{dataconn}->get_collection( 'biosamples' );
  $cursor	      =		$bscoll->find( { id => { '$in' => [ keys %biosIds ] } } )->fields( { id => 1, $groupAttr => 1 } );
  my $bioS      =   [ $cursor->all ];

  for my $i (0..$#{ $pgx->{samples} }) {

    my $csId    =   $pgx->{samples}->[$i]->{id};
    my $bsId    =   $pgx->{samples}->[$i]->{biosample_id};
    
    $pgx->{samples}->[$i]->{sortkey}    =  'NA';
    $pgx->{samples}->[$i]->{sortlabel}  =  'not specified';
   
    if ($bsId !~ /.../) { next }    
    my ($thisBios)  =   grep{ $_->{id} eq $bsId } @$bioS;
    my ($thisbioc)  =   grep{ $_->{type}->{id} =~ /$groupType/ } @{ $thisBios->{$groupAttr} };
    if ($thisbioc->{type}->{id} !~ /.../) { next }
        
    $pgx->{samples}->[$i]->{sortkey}    =   $thisbioc->{type}->{id};
    $pgx->{samples}->[$i]->{sortlabel}  =   ( $thisbioc->{type}->{label} =~ /.../ ? $thisbioc->{type}->{label} : $thisbioc->{type}->{id} );  
    $pgx->{samples}->[$i]->{sortlabel}  =~  s/^.+?\:\:?//;  

  }  

  return $pgx;
  
}

################################################################################

sub pgx_callset_labels_from_file {

  my $pgx       =   shift;
  my $sortFile  =   shift;
  
  my $fallbackK =   'NA';
  my $fallbackL =   'not specified';

  if (! -f $sortFile)  { return $pgx }

  # this assumes that the first column contains an entry for the selected id (or id)
  # the second then the associated label
  my $customSort  =   {};

  my $table     =   read_file_to_split_array($sortFile);
    foreach (@$table) {
      if ($_->[0] =~ /\w\w/ && $_->[1] =~ /\w/) {
        $customSort->{$_->[0]}  =   {
          sortkey     =>  $_->[1] =~ /\w\w/ ? $_->[1] : $fallbackK,
          sortlabel   =>  $_->[2] =~ /\w\w/ ? $_->[2] : $fallbackL,
        };
  }}
  
  for my $i (0..$#{ $pgx->{samples} }) {

    my $csId    =   $pgx->{samples}->[$i]->{id};
    
    $pgx->{samples}->[$i]->{sortkey}    =   $fallbackK;
    $pgx->{samples}->[$i]->{sortlabel}  =   $fallbackL;
    
    if ($customSort->{$csId}->{sortkey} !~ /.../) { next }    
        
    $pgx->{samples}->[$i]->{sortkey}    =   $customSort->{$csId}->{sortkey};
    $pgx->{samples}->[$i]->{sortlabel}  =   $customSort->{$csId}->{sortlabel};
    
  }

  return $pgx;
  
}

1;
