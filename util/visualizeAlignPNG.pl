#!/usr/local/bin/perl -w
use strict;

#
# Module Dependence
#
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp qw/ tempfile tempdir /;
use MultAln;
use GD::Simple;

#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use CrossmatchSearchEngine;
use SearchResult;
use SearchResultCollection;

my $config           = $RepModelConfig::configuration;
my $ucscToolsDir =  $config->{'UCSCTOOLS_DIR'}->{'value'};

my @getopt_args = (
                    '-genome=s',
                    '-outfile=s',
                    '-fasta=s',
                    '-coreoutfile=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  die;
}

my $genome_file;

if ( -S "genome" ) {
  $genome_file = "genome";
}
if ( exists $options{'genome'} )
{
  $genome_file = $options{'genome'};
}else {
  warn "Genome not provided.  No flanking data will be displayed\n";
}

my $outFile;
if ( -S "out" ) {
  $outFile = "out";
}
if ( exists $options{'outfile'} ) {
  $outFile = $options{'outfile'};
}

# Obtain seq lengths -- do this once
my %seqLens = ();
if ( $genome_file ) 
{
  open IN, "$ucscToolsDir/twoBitInfo $genome_file stdout|"
    or die "Could not run $ucscToolsDir/twoBitInfo on $genome_file!\n";
  my $lines = 0;
  while ( <IN> )
  {
    $lines++;
    if ( /^(\S+)\s+(\d+)/ )
    {
      $seqLens{$1} = $2;
    }
  }
  close IN;
}

my $maxLen = 0;
my $maxSeqID;
my $maxStart;
my $maxEnd;
my $mAlign;
if ( exists $options{'outfile'} ) {
  my $resultCollection =
    CrossmatchSearchEngine::parseOutput( searchOutput => $outFile );

  if ( $options{'coreoutfile'} ) {
    my $coreResultCollection =
      CrossmatchSearchEngine::parseOutput( searchOutput => $options{'coreoutfile'} );
    for ( my $i = 0; $i < $coreResultCollection->size(); $i++ ){
       my $result = $coreResultCollection->get($i); 
       if ( $result->getQueryEnd() - $result->getQueryStart() + 1 > $maxLen ) {
         $maxSeqID = $result->getQueryname();
         $maxStart = $result->getQueryStart();
         $maxEnd =  $result->getQueryEnd();
         $maxLen = $result->getQueryEnd() - $result->getQueryStart() + 1;
       }
    }
    undef $coreResultCollection;
  }
  
  $mAlign = MultAln->new(
                             referenceSeq              => "",
                             searchCollection          => $resultCollection,
                             searchCollectionReference => MultAln::Subject
                        );
}elsif ( exists $options{'fasta'} ) {
  my @seqs;
  my $seq;
  my $id;
  open my $IN, "<$options{'fasta'}" or die "Could not open $options{'fasta'} for reading";
  # Simple FASTA reader
  my %idHash = ();
  while (<$IN>) {
    if ( /^>(\S+)/ )
    {
      my $tmpID = $1;
      if ( defined $idHash{$tmpID} ) {
        my $ver = 1;
        while ( defined $idHash{$tmpID . "_$ver"} )
        {
          $ver++;
        }
        warn "WARN File contains a duplicate identifier \"$tmpID\".  A suffix of \"_$ver\"\n" .
             "     will be appended to this occurence.\n";
        $tmpID = $tmpID . "_$ver";
      }
      $idHash{$tmpID}++;
      if ( $seq )
      {
        $seq = uc($seq);
        # Convert prefix/suffix "-"s to spacers
        #if ( $seq =~ /^(\-+)/ ){
        #  substr($seq,0,length($1)) = "."x(length($1));
        #}
        #if ( $seq =~ /(\-+)$/ ) {
        #  substr($seq,length($seq)-length($1)-1) = "."x(length($1));
        #}
        $seq =~ s/\-/./g;
        push @seqs, [ " ", $seq, " " ];
      }
      $seq = "";
      $id = $tmpID;
      next;
    }
    s/[\s\n\r]+//g;
    $seq .= $_;
  }
  if ( $seq )
  {
    # Convert prefix/suffix "-"s to spacers
    #if ( $seq =~ /^(\-+)/ ){
    #  substr($seq,0,length($1)) = "."x(length($1));
    #}
    #if ( $seq =~ /(\-+)$/ ) {
    #  substr($seq,length($seq)-length($1)-1) = "."x(length($1));
    #}
    $seq =~ s/\-/./g;
    $seq = uc($seq);
    push @seqs, [ " ", $seq, " " ];
  }
  close $IN;
  getMultAlignPNG( seqs => \@seqs );
  exit;
  ###
}else {
  die "Must provide either -outfile or -fasta";
}
  
if ($genome_file eq "" ) {
  die "At this time genome is not an optional argument when used with -outfile";
}
my $genome     = $genome_file;
my $flankleft  = 50;
my $flankright = 50;
my $DEBUG      = 0;

my @ranges = ();
my @order  = ();
my ( $tmpFH, $tmpFilename ) =
    tempfile( UNLINK => 0, SUFFIX => ".seqlist", DIR => "." );
my $maxAlignedLen = 0;
for ( my $i = 0 ; $i < $mAlign->getNumAlignedSeqs() ; $i++ )
{
  my $seq         = $mAlign->getAlignedSeq( $i );
  my $raw_id      = $mAlign->getAlignedName( $i );
  my $align_start = $mAlign->getAlignedSeqStart( $i );
  my $align_end   = $mAlign->getAlignedSeqEnd( $i );
  my $orient      = $mAlign->getAlignedOrientation( $i );
  my $seqOffset   = $mAlign->getAlignedStart( $i );
  if ( length( $seq ) + $seqOffset > $maxAlignedLen ) {
    $maxAlignedLen = length( $seq ) + $seqOffset;
  }

  #print "$raw_id $align_start,$align_end:$orient : $seq\n";

  my $align_remain = 0;

  my $id = "";
  my $start;
  my $end;

  # These have traditionally been in 1-based full-closed coordinates
  if ( $raw_id =~ /^(\S+)\_(\d+)\_(\d+)\_?([R\+\-]?)$/ )
  {
    $id    = $1;
    $start = $2;
    $end   = $3;
  } elsif ( $raw_id =~ /^(\S+):(\d+)-(\d+)\_?([R\+\-]?)$/ )
  {
    $id    = $1;
    $start = $2;
    $end   = $3;
  } else
  {
    die "I don't know how to parse this id: $raw_id\n";
  }
  $ranges[ $i ] = { 'id' => $id, 'seq' => $seq };

  if ( !exists $seqLens{$id} )
  {
    warn "Could not find $id in 2bit file!\n";
  }
  my $slen = $seqLens{$id};

  # Does the current sequence orientation still make sense or has
  # it reverted back to the forward strand?
  if (    $raw_id =~ /_R/ || $raw_id =~ /_-$/ 
       || $orient eq "-" && !( $raw_id =~ /_R/ && $orient eq "-" ) )
  {

    # Reversed orientation ( overall )
    # 1-based, fully closed
    $ranges[ $i ]->{'orient'} = "-";
    my $leftStart = $end - $align_start;
    my $leftEnd   = $leftStart - 1;
    if ( $leftEnd + $flankleft > $slen )
    {
      $leftEnd = $slen;
    } else
    {
      $leftEnd += $flankleft;
    }
    if ( $leftEnd - $leftStart > 0 )
    {
      $ranges[ $i ]->{'leftStart'} = $leftStart;
      $ranges[ $i ]->{'leftEnd'}   = $leftEnd;
      push @order, [ "L", $i ];

      #push @ranges, [ $i, $id, "L", $leftStart, $leftEnd, "-" ];
      print $tmpFH "$id:" . ( $leftStart - 1 ) . "-" . $leftEnd . "\n";
    }

    my $rightEnd   = $end - $align_end + 2;
    my $rightStart = $rightEnd + 1;
    if ( $rightStart - $flankright < 1 )
    {
      $rightStart = 1;
    } else
    {
      $rightStart -= $flankright;
    }
    if ( $rightEnd - $rightStart > 0 )
    {
      $ranges[ $i ]->{'rightStart'} = $rightStart;
      $ranges[ $i ]->{'rightEnd'}   = $rightEnd;
      push @order, [ "R", $i ];

      #push @ranges, [ $i, $id, "R", $rightStart, $rightEnd, "-" ];
      print $tmpFH "$id:" . ( $rightStart - 1 ) . "-" . $rightEnd . "\n";
    }

    #print "Extension: $leftStart-$leftEnd, ... ,$rightStart-$rightEnd (-)\n";
  } else
  {

    # Normal left/right
    $ranges[ $i ]->{'orient'} = "+";
    my $leftEnd = $start + $align_start -
        2;    # -2 = shift to left of alignment and correct for 1-based math.
    my $leftStart = $leftEnd + 1;
    if ( $leftStart - $flankleft < 1 )
    {
      $leftStart = 1;
    } else
    {
      $leftStart -= $flankleft;
    }
    if ( $leftEnd - $leftStart > 0 )
    {
      $ranges[ $i ]->{'leftStart'} = $leftStart;
      $ranges[ $i ]->{'leftEnd'}   = $leftEnd;
      push @order, [ "L", $i ];

      #push @ranges, [ $i, $id, "L", $leftStart, $leftEnd, "+" ];
      print $tmpFH "$id:" . ( $leftStart - 1 ) . "-" . $leftEnd . "\n";
    }

    my $rightStart = $end + $align_start;
    my $rightEnd   = $rightStart - 1;
    if ( $rightEnd + $flankright > $slen )
    {
      $rightEnd = $slen;
    } else
    {
      $rightEnd += $flankright;
    }
    if ( $rightEnd - $rightStart > 0 )
    {
      $ranges[ $i ]->{'rightStart'} = $rightStart;
      $ranges[ $i ]->{'rightEnd'}   = $rightEnd;
      push @order, [ "R", $i ];

      #push @ranges, [ $i, $id, "R", $rightStart, $rightEnd, "+" ];
      print $tmpFH "$id:" . ( $rightStart - 1 ) . "-" . $rightEnd . "\n";
    }

    #print "Extension: $leftStart-$leftEnd, ... ,$rightStart-$rightEnd (+)\n";
  }
}
close $tmpFH;

# Generate sequences
my $cmd = "$ucscToolsDir/twoBitToFa -seqList=$tmpFilename $genome stdout";
open IN, "$cmd|" or die "Could not run $cmd!\n";
my $idx = 0;
my $seq = "";
my $id;
my $start;
my $end;

while ( <IN> )
{
  if (    />(\S+)\:(\d+)-(\d+)/
       || />(\S+)/ )
  {
    my $tmp_id    = $1;
    my $tmp_start = $2;
    my $tmp_end   = $3;
    $tmp_start = 0 if ( !defined $tmp_start );
    $tmp_end = $seqLens{$tmp_id} if ( !defined $tmp_end );
    if ( $seq )
    {
      my $leftOrRight = $order[ $idx ][ 0 ];
      my $rangeIdx    = $order[ $idx ][ 1 ];
      my $outID       = $id . "_" . ( $start + 1 ) . "_$end";
      if ( $ranges[ $rangeIdx ]->{'orient'} eq "-" )
      {
        $seq = reverse( $seq );
        $seq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;
        $outID .= "_R";
      }
      if ( $leftOrRight eq "R" )
      {
        $ranges[ $rangeIdx ]->{'rightFlank'} = $seq;
      } else
      {
        $ranges[ $rangeIdx ]->{'leftFlank'} = $seq;
      }
      $idx++;
    }
    $id    = $tmp_id;
    $start = $tmp_start;
    $end   = $tmp_end;
    $seq   = "";
    next;
  }
  s/[\n\r\s]+//g;
  $seq .= uc( $_ );
}
if ( $seq )
{
  my $leftOrRight = $order[ $idx ][ 0 ];
  my $rangeIdx    = $order[ $idx ][ 1 ];
  my $outID       = $id . "_" . ( $start + 1 ) . "_$end";
  if ( $ranges[ $rangeIdx ]->{'orient'} eq "-" )
  {
    $seq = reverse( $seq );
    $seq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;
    $outID .= "_R";
  }
  if ( $leftOrRight eq "R" )
  {
    $ranges[ $rangeIdx ]->{'rightFlank'} = $seq;
  } else
  {
    $ranges[ $rangeIdx ]->{'leftFlank'} = $seq;
  }
}
close IN;

#print "" . Dumper( \@ranges ) . "\n";
#exit;

my @seqs = ();
for ( my $i = 0 ; $i < $mAlign->getNumAlignedSeqs() ; $i++ )
{
  my $seq         = $mAlign->getAlignedSeq( $i );
  my $raw_id      = $mAlign->getAlignedName( $i );
  my $align_start = $mAlign->getAlignedSeqStart( $i );
  my $align_end   = $mAlign->getAlignedSeqEnd( $i );
  my $orient      = $mAlign->getAlignedOrientation( $i );
  my $seqOffset   = $mAlign->getAlignedStart( $i );

  if ( $raw_id eq $maxSeqID ) {
    print "FOUND MaxID\n";
    if ( $align_start <= $maxStart && $align_end >= $maxEnd ) {
      print "And it's within range!\n";
    }
  }

  my $leftFlank = "";
  if ( exists $ranges[ $i ]->{'leftFlank'} )
  {
    $leftFlank = $ranges[ $i ]->{'leftFlank'};
  }
  $leftFlank = " " x ( $flankleft - length( $leftFlank ) ) . $leftFlank;
  my $rightFlank = "";
  if ( exists $ranges[ $i ]->{'rightFlank'} )
  {
    $rightFlank = $ranges[ $i ]->{'rightFlank'};
  }
  $rightFlank = $rightFlank . " " x ( $flankright - length( $rightFlank ) );

  $seq = " " x ( $seqOffset ) . $seq;
#print "i=$i length(seq)=" . (length($seq))." maxAlignedLen=$maxAlignedLen\n";
  $seq = $seq . " " x ( $maxAlignedLen - length( $seq ) );

  push @seqs, [ $leftFlank, $seq, $rightFlank ];

  #  print "$raw_id $align_start,$align_end:$orient : $seq\n";
}

#@seqs = sort { $a->[2] cmp $b->[2] } @seqs;
@seqs = sort { reverse( $a->[ 0 ] ) cmp reverse( $b->[ 0 ] ) } @seqs;

#print "" . Dumper( \@seqs ) . "\n";
#exit;

getMultAlignPNG( seqs => \@seqs );

unless ( $DEBUG )
{
  unlink( $tmpFilename ) if ( -e $tmpFilename );
}

########################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##

=head2 getMultAlignPNG()

  Use:  my $svg = getMultAlignPNG( );

  Returns a PNG object containing the ...

=cut

##-------------------------------------------------------------------------##
sub getMultAlignPNG
{
  my %namedArgs = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  die "$subroutine: Missing \"seqs\" argument!\n"
      unless ( defined $namedArgs{'seqs'} );
  my @seqs = @{ $namedArgs{'seqs'} };

  # For DNA only right now
  # AliView
  my %baseColors = (
                     'T' => "LightSkyBlue",
                     'C' => "LightGreen",
                     'G' => "Gold",
                     'A' => "IndianRed",
                     '.' => "white",
                     '-' => "white",
                     ' ' => "white"
  );


 # UniPro
 #my %baseColors = ( 'T' => "red", 'C' => "lime", 'G' => "blue", 'A' => "gold");

  #my @seqs = (
  #             "GATAGCTATGAC"
  #                 . "TTA..CTT.AGCTAGGATTAG.ATCA...GTGCAGGTA"
  #                 . "AGTCAGACTAGACGTTCAG",
  #             "TAGGAATGCTGC"
  #                 . "TAAGGCCC.AGCTAAGATTAG.ATCA...GTGCTTGTA"
  #                 . "TTAGGCTACGATGCGTGAT",
  #             "ACAGGGGGATAT"
  #                 . "TTA..CCC.AGCTAGGATTAC.ATCA...GTGCAGGTA"
  #                 . "GATAATCGGATACTATCGC",
  #             "TTCATACCCCCT"
  #                 . "TTAG.CGC.AGGTAGGATTAG.ATCAAAAGAGCAGGTA"
  #                 . "AATACGTGTCGAGCATGTT",
  #             "ACCATTTGGGTA"
  #                 . "TTAG.CCC.AGCTAGGTTTAG.ATCA...GTGCAGGCA"
  #                 . "TTAGCAGGAGTACGCGTTA",
  #             "CGTAGGTATGGC"
  #                 . "TTAGGCCCAAGCTAGGATTACCATCA...GTGAAGGTA"
  #                 . "ACTAGTTGGGGATTAGAAA"
  #);

  #
  # Data Characteristics
  #
  my $numAlignments = $#seqs + 1;
  my $flankLength = length($seqs[0][0]);
  my $alignLength   = 0;
  for ( my $i = 0 ; $i <= $#seqs ; $i++ )
  {
    my $combLen =
        length( $seqs[ $i ][ 0 ] ) + length( $seqs[ $i ][ 1 ] ) +
        length( $seqs[ $i ][ 2 ] );
    $alignLength = $combLen if ( $combLen > $alignLength );
  }


  #
  # Visual constants and relativities
  #
  my $vizMultAlignFontSize = 8;

  my $topMargin    = 10;
  my $leftMargin   = 10;
  my $rightMargin  = 10;
  my $bottomMargin = 10;
  my $sectionMargin = 50;

  my $alignedPosHeight = 5;
  my $alignedPosWidth  = 5;
  my $heatMapWidth = $alignLength * $alignedPosWidth;
  my $heatMapHeight = $numAlignments * $alignedPosHeight;

  my $textHeight = 10;
  my $textWidth = 8;
  my $textAlignWidth = 2*((10*$textWidth) + 20 + (10*$textWidth));
  my $textAlignHeight = $numAlignments * $textHeight;
  my $totalPNGWidth  = sprintf( "%0.0f", $leftMargin + $heatMapWidth + $rightMargin );
  my $totalPNGHeight = sprintf( "%0.0f", $topMargin + $heatMapHeight + $sectionMargin + $textAlignHeight + $bottomMargin );

  # Set origin to be inside the axis
  my $graphOriginX = $leftMargin;
  my $graphOriginY = $topMargin;

  my $img = GD::Simple->new( $totalPNGWidth + 1000, $totalPNGHeight + 100 );
  my $black = $img->colorAllocate(0,0,0);
  my %fontColors = (
                     'T' => $img->colorAllocate(0,0,255),
                     'C' => $img->colorAllocate(0,255,0),
                     'G' => $img->colorAllocate(255,155,48),
                     'A' => $img->colorAllocate(255,0,0),
                     'N' => $img->colorAllocate(0,0,0),
                     '.' => $img->colorAllocate(0,0,0),
                     '-' => $img->colorAllocate(0,0,0),
                     ' ' => $img->colorAllocate(0,0,0)
  );

  my $boundaryTitleBarHeight = 36;
  $img->string(gdGiantFont, $leftMargin, $topMargin + $heatMapHeight + $sectionMargin, "Left Alignment Boundary", $black);
  $img->string(gdGiantFont, $leftMargin, $topMargin + $heatMapHeight + $sectionMargin + 18, "Flanking", $black);
  $img->string(gdGiantFont, $leftMargin + $sectionMargin + (10*12), $topMargin + $heatMapHeight + $sectionMargin + 18, "Aligned", $black);

  $img->string(gdGiantFont, $leftMargin + (3*$sectionMargin) + (20*12) , $topMargin + $heatMapHeight + $sectionMargin, "Right Alignment Boundary", $black);
  $img->string(gdGiantFont, $leftMargin + (3*$sectionMargin) + (20*12), $topMargin + $heatMapHeight + $sectionMargin + 18, "Flanking", $black);
  $img->string(gdGiantFont, $leftMargin + (6*$sectionMargin) + (20*12),  $topMargin + $heatMapHeight + $sectionMargin + 18, "Aligned", $black);

  my $y = $graphOriginY;
  my $textY = $topMargin + $heatMapHeight + $sectionMargin + $boundaryTitleBarHeight;
  for ( my $i = 0 ; $i < $numAlignments ; $i++ )
  {
    my $x = $graphOriginX;
    for ( my $j = 0 ; $j < $alignLength ; $j++ )
    {
      my $combSeq = $seqs[ $i ][ 0 ] . $seqs[ $i ][ 1 ] . $seqs[ $i ][ 2 ];
      my $base    = " ";
      if ( $j < length( $combSeq ) )
      {
        $base = substr( $combSeq, $j, 1 );
      }
      my $color = "white";
      $color = $baseColors{$base} if ( exists $baseColors{$base} );
      $img->bgcolor( $color );
      $img->fgcolor( $color );
      $img->rectangle( $x, $y, $x + $alignedPosWidth, $y + $alignedPosHeight );
      
      $x += $alignedPosWidth;
    }
    $y += $alignedPosHeight;

    my $textX = $leftMargin;
    my $skipped = 0;
    if ( substr($seqs[$i][1],0,10) =~ /[ACGTN]/ ){
    # Left Flanking
    for ( my $k = 10; $k > 0; $k-- )
    {
      my $base = substr($seqs[$i][0], -$k, 1);
      $img->string(gdLargeFont, $textX , $textY, $base, $fontColors{$base});
      $textX += 12;
    }
    $textX += $sectionMargin;
    # Left Aligned
    for ( my $k = 0; $k < 10; $k++ )
    {
      my $base = substr($seqs[$i][1], $k, 1);
      $img->string(gdLargeFont, $textX, $textY, $base, $fontColors{$base});
      $textX += 12;
    }
    }else { $textX += 240+$sectionMargin; $skipped++; }

    $textX += 2*$sectionMargin;

    if ( substr($seqs[$i][1],-10,10) =~ /[ACGTN]/ ) {
    # Right Aligned
    for ( my $k = 10; $k > 0; $k-- )
    {
      my $base = substr($seqs[$i][1], -$k, 1);
      $img->string(gdLargeFont, $textX , $textY, $base, $fontColors{$base});
      $textX += 12;
    }
    $textX += $sectionMargin;
    # Right Flanking
    for ( my $k = 0; $k < 10; $k++ )
    {
      my $base = substr($seqs[$i][2], $k, 1);
      $img->string(gdLargeFont, $textX, $textY, $base, $fontColors{$base});
      $textX += 12;
    }
    }else { $skipped++; }

    if ( $skipped < 2 ) {
    $textY += 18;
    }

  }
  # Left flank divider
  $img->line( $leftMargin + ($flankLength * $alignedPosWidth), $topMargin, 
              $leftMargin + ($flankLength * $alignedPosWidth), $topMargin + ($numAlignments*$alignedPosHeight), $black );
  # Right flank divider
  $img->line( $leftMargin + ($flankLength * $alignedPosWidth) + ($maxAlignedLen * $alignedPosWidth), $topMargin, 
              $leftMargin + ($flankLength * $alignedPosWidth) + ($maxAlignedLen * $alignedPosWidth), $topMargin + ($numAlignments*$alignedPosHeight), $black );


  open my $out, '>', 'img.png' or die;
  binmode $out;
  print $out $img->png;

  return;

}

1;
