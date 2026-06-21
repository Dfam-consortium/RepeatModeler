#!/usr/bin/env perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) resolveIndels.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Arian Smit asmit@systemsbiology.org.
##  Description:
##      A revised Blocker/AutoRunBlocker script to resolve indels in
##      tranitively aligned MSAs.
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
###############################################################################
#
#

=head1 NAME

resolveIndels - Remove rare indels in the MSA consensus

=head1 SYNOPSIS

  resolveIndels [-version] 
                 [-ruzzo_tompa_threshold #]
                 ( [-discrete_windows <#,#,#,...>] |
                   [-window_min #] [-window_max #] )
                 [-aggregation_method [tile|cluster]]
                 [-min_gap #]
                 [-min_copy #]
                 [-min_ratio #]
                 [-cons <output consensus filename>]
                 [-interactive]
                 -msa <Linup *.ali|*.fasta|*.multaln|*.out>

=head1 DESCRIPTION

A jacknife tool for identifying rare indels in transitively
aligned MSA.  


The first step in this process is to identify probable indel
column ranges in the alignment.  There are two methods 
available here.  The first is the Blocker approach - run 
a window across the MSA and evaluate the sequence length
distribution in each window.  The window may be a single
discrete size, a set of discrete sizes, or a range:

To run with discrete window sizes 7,15,24,5:

./resolveIndels -discrete_windows 7,15,24,5 -msa foo.ali

To run with windows from 5-24bp:

./resolveIndels -window_min 5 -window_max 24 -msa foo.ali

Alternatively the windows are obtained in using the same
algorithm used by Refiner (Ruzzo-Tompa find all low scoring
sub subalignments):

./resolveIndels -ruzzo_tompa_threshold 1 -msa foo.ali

Once the candidate ranges have been identified ( reported as
"Initial Ranges" by the tool ), the ranges are then either
tiled or clustered.  Here we generate a tiling path requiring
at least 10bp between ranges:

./resolveIndels -min_gap 10 -aggregation_method tile 
                -discrete_windows 7,15,24,5 -msa foo.ali

Alternatively, the ranges can be clustered using:

./resolveIndels -aggregation_method cluster
                -discrete_windows 7,15,24,5 -msa foo.ali


The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2012-2021 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>
Arian Smit <asmit@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use Getopt::Long;
use Data::Dumper;
use Cwd;
use File::Spec;
use File::Basename;
use Carp;

# RepeatModeler Libraries
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/..";
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use MultAln;
use SeedAlignment;

# RepeatMasker Libraries
use RepeatUtil;
use SearchResult;
use SearchResultCollection;
use WUBlastSearchEngine;
use NCBIBlastSearchEngine;
use CrossmatchSearchEngine;
use FastaDB;

my $Version    = $RepModelConfig::VERSION;
my $ucscToolsDir = "/usr/local/bin";

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;
$DEBUG = 1 if ( $RepModelConfig::DEBUG == 1 );

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',          # print out the version and exit
                    '-msa=s',
                    '-aggregation_method=s',
                    '-ruzzo_tompa_threshold=s',
                    '-discrete_windows=s',
                    '-interactive',
                    '-window_min=s',
                    '-window_max=s',
                    '-min_gap=s',
                    '-min_copy=s',
                    '-min_ratio=s',
                    '-cons=s',
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) )
{
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0 | less";
  exit( 1 );
}


if ( $options{'version'} )
{
  print "$Version\n";
  exit;
}

my $aggregationMethod = "tile";
if ( $options{'aggregation_method'} ) {
  if ( $options{'aggregation_method'} eq "cluster" ) {
    $aggregationMethod = "cluster";
  }elsif ( $options{'aggregation_method'} eq "tile" ) {
    # already default
  }else {
    print "Unknow aggregation method $options{'aggregation_method'}, please use 'tile' or 'cluster'.\n";
  }
}
 
my $minSep = 0;
$minSep = $options{'min_gap'} if ( $options{'min_gap'} );

my $method = "blocker";
if ( $options{'discrete_windows'} || $options{'window_min'} || $options{'window_max'} ) {
  $method = "blocker";
}elsif ( $options{'ruzzo_tompa_threshold'} )  {
  $method = "ruzzo_tompa";
}

if ( $options{'discrete_windows'} && ( $options{'window_min'} || $options{'window_max'} )) {
  print "\n\nOption discrete_windows cannot be combined with window_min/window_max\n\n";
  usage()
}

if ( ! ( $options{'window_min'} && $options{'window_max'} ) &&
     ( $options{'window_min'} || $options{'window_max'} ) ) {
  print "\n\nOptions window_min and window_max must be used together.\n\n";
  usage()
}

## TODO: Make this a parameter
open (IN, "$FindBin::RealBin/../Matrices/linupmatrix ") or die;
my @bases;
my @l;
my %matrix;
while (<IN>) {
  s/^\s+//;
  if (/^A/) {
    @l = split;
    @bases = @l;
  } elsif (/^[\d-]/) {
    my @s = split;
    my $l = shift @bases;
    for (my $i = 0; $i <= $#s; ++$i) {
      $matrix{$l}{$l[$i]} = $s[$i];
    }
  }
}
close IN;
 
my @windowSizes = ();
my $methodLabel = "";
if ( $options{'discrete_windows'} ) {
  @windowSizes = split(/\s*,\s*/,$options{'discrete_windows'});
  $methodLabel = "Discrete window sizes " . join(",",@windowSizes);
}elsif ( $options{'window_min'} && $options{'window_max'} ) {
  for ( my $i = $options{'window_min'}; $i <= $options{'window_max'}; $i++ ) {
    push @windowSizes, $i;
  }
  $methodLabel = "Continuous window sizes " . $options{'window_min'} . " to " . $options{'window_max'};
}elsif ( $method eq "ruzzo_tompa" ) {
  $methodLabel = "Ruzzo Tompa Threshold : " . $options{'ruzzo_tompa_threshold'} . "\n";
}else {
  print "\n\nMust supply either option discrete_windows or both window_min and window_max or ruzzo_tompa_threshold\n\n";
  usage();
}

my $inputFile = $options{'msa'};
if ( !-s $inputFile )
{
  print "\nCannot locate file!: $inputFile\n\n";
  usage();
}

my $inputFileDir = cwd;
my($filename, $dirs, $suffix) = fileparse($inputFile);
$inputFileDir = $dirs;

my $copymin = 4;
$copymin = $options{'min_copy'} if ( $options{'min_copy'} );
my $ratio = 2;
$ratio = $options{'min_ratio'} if ( $options{'min_ratio'} );
my $verbose = 1;

my ($mAlign, $fileType) = RepeatUtil::openAsMultAln( $inputFile );
my $cons = $mAlign->consensus();

# Create an array to determine consensus position from msa column
# in 1-based coordinates.
#   ie.
#              cons: A--GTTA-TTT--
#      msaToConsPos: 1112345567888
#
my @msaToConPos = ();
my $cIdx = 0;
for ( my $i=0; $i < length($cons); $i++ ) {
  my $char = substr($cons,$i,1);
  $cIdx++ if ( $char ne "-" );
  push @msaToConPos, $cIdx;
}

if ( $verbose ) {
  print STDERR "##\n";
  print STDERR "## resolveIndels.pl\n";
  print STDERR "##\n";
  print STDERR "# MSA File = $inputFile\n";
  print STDERR "#  - MSA Columns = " . length($cons) . "\n";
  my $gaplessCons = $cons;
  $gaplessCons =~ s/-//g;
  print STDERR "#  - MSA Consensus Length = " . length($gaplessCons) . "\n";
  print STDERR "# Min Sequences with Alternative Length (min_copy) = $copymin\n";
  print STDERR "# Min Ratio of Alternative Length to Current Length Blocks (min_ratio) = $ratio\n";
  print STDERR "# Aggregation Method = $aggregationMethod\n";
  print STDERR "# Min Gap Between Ranges = $minSep\n";
  print STDERR "#\n";
}


## Use ruzzo_tompa
my $ranges;
if ( $method eq "ruzzo_tompa" ) {
  my $threshold = $options{'ruzzo_tompa_threshold'};
  if ( $verbose ) {
    print STDERR "# Inverted Ruzzo Tompa Method\n";
    print STDERR "# ---------------------------\n";
    print STDERR "# Threshold = $threshold\n";
    print STDERR "# Scoring Matrix = comparison.matrix (symmetric, avg +10 diagonal)\n";
    print STDERR "# GapOpen/GapExt = -40/-15\n";
    print STDERR "# NOTE: These are low scoring sub-alignment ranges.  They are not necessarily characterized by\n";
    print STDERR "#       length variations.  Therefore, only the ones with length variation passing the ratio\n";
    print STDERR "#       set by the program (marked with '*' will be passed on for cluster/tiling\n";
  }
  print "Ruzzo-Tompa Ranges (by MSA column indices)\n";
  print "  Start\tEnd\tScore\n";
  print "-------------------\n";
  my ($ruzzoTompaBlocks,$ruzzoTompaScores) = $mAlign->getLowScoringAlignmentColumns( 'threshold' => $threshold );

  # DEBUG:
  if ( 0 ) {
    my $html = getMultAlnHTML(
      multAln => $mAlign,
      headersAboveMSA => [
         ['lowQualScore', $ruzzoTompaScores]
      ],
      coloredBlocks => [
        ['#CD5C5C', $ruzzoTompaBlocks],
      ],
      alignmentBlocks => $ruzzoTompaBlocks,
      inclRef => 1,
      DEBUG => 1,
    );
    open OUT,">/home/rhubley/public_html/ri.html" or die;
    print OUT $html;
    close OUT;
  }

  $ranges = [];
  for ( my $i = 0; $i <= $#{$ruzzoTompaBlocks}; $i++ ) {
    my $block = $ruzzoTompaBlocks->[$i];
    my $score = $ruzzoTompaScores->[$block->[0]];
    my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
             &evalMSABlock( $mAlign, $block->[0], $block->[1]-1, $copymin, $ratio, \%matrix, $cons);
    if ( $numMajorityLen ) {
      push @{$ranges}, [$score, $block->[0], $block->[1]-1];
      print "  " . $block->[0] . "\t" . ($block->[1]-1) . "\t$score\t*\n";
    }else { 
      print "  " . $block->[0] . "\t" . ($block->[1]-1) . "\t$score\n";
    }
  }
}elsif ( $method eq "blocker" ) {
  if ( $verbose ) {
    print STDERR "# Extended Blocker Method\n";
    print STDERR "# ---------------------------\n";
    print STDERR "# $methodLabel\n";
    print STDERR "#\n";
  }
  # Grab all windows as specified by user
  $ranges = &identifyUnstableBlockerWindows( $mAlign, \@windowSizes, $copymin, $ratio, \%matrix );


}

# Sort by ratio position
my @sortedRanges = sort { $b->[0] <=> $a->[0] } @{$ranges};

# Just print them out
print "Initial Ranges:\n";
print "  Ratio   Range      Consensus\n";
print "  -------------------------------\n";
foreach my $range ( @sortedRanges ) {
  print "" . sprintf("%7.2f",$range->[0]) . " : " . $range->[1] . " - " . $range->[2] . ", " . $range->[3] . "\n";
}
print "\n";

# Create a tiling path of high scoring ranges seperated by $minSep bases
my @tilingPath = ();

if ( $aggregationMethod eq "tile" ) {
  if ( $minSep > 0 ) {
    print "Tiling Path: min separation $minSep bp\n";
  }elsif ( $minSep == 0 ) {
    print "Tiling Path: non-overlaping\n";
  }else {
    print "Tiling Path: max overlap" . abs($minSep) . "\n";
  }
  foreach my $range ( @sortedRanges ) {
    my $overlap = 0;
    my $rangeID = 0;
    foreach my $existingRange ( @tilingPath ) {
      $rangeID++;
      if ( ( $range->[1] >= $existingRange->[1]-$minSep && $range->[1] <= $existingRange->[2]+$minSep ) ||
           ( $range->[2] >= $existingRange->[1]-$minSep && $range->[2] <= $existingRange->[2]+$minSep ) ||
           ( $range->[1] < $existingRange->[1]-$minSep && $range->[2] > $existingRange->[2]+$minSep) ) {
        # Overlaps disqualify
        $overlap = $rangeID;
        last;
      }
    }
    unless ( $overlap ) {
      push @tilingPath, $range;
      push @{$range}, "*";
      push @{$range}, scalar(@tilingPath);
    }else {
      push @{$range}, " ";
      push @{$range}, $overlap;
    }
  }
  print "  Sel? Ratio  MSA_Range  Cons_Range Delta       Consensus\n"; 
  print "  ---------------------------------------------------------------------\n";
  @sortedRanges = sort { $a->[1] <=> $b->[1] || $b->[2] <=> $a->[2] } @{$ranges};
  my $grpIdx = 0;
  my $minRange;
  my $maxRange;
  my @colWidths;
  for ( my $i = 0; $i <= $#sortedRanges; $i++ ) {
    my $range = $sortedRanges[$i];
    
    # Segment tilling group and look-ahead to get min/max msa col range
    if ( $range->[5] != $grpIdx ) {
      $minRange = 99999999;
      $maxRange = -99999999;
      $grpIdx = $range->[5];
      for ( my $j = $i; $j <= $#sortedRanges; $j++ ) {
        my $srange = $sortedRanges[$j];
        last if ( $srange->[5] != $grpIdx );
        $minRange = $srange->[1] if ( $srange->[1] < $minRange );
        $maxRange = $srange->[2] if ( $srange->[2] > $maxRange );
      }
      my $consSegment = substr($cons,$minRange, $maxRange-$minRange+1);
      $consSegment =~ s/-//g;
      $colWidths[0] = length($minRange);
      $colWidths[1] = length($maxRange);
      $colWidths[2] = length($msaToConPos[$minRange]);
      $colWidths[3] = length($msaToConPos[$maxRange]);
      print "  Region Consensus:" . " "x($colWidths[0] + $colWidths[1]) . 
            sprintf("%$colWidths[2]d",$msaToConPos[$minRange]) .  " - " . sprintf("%$colWidths[3]d",$msaToConPos[$maxRange]) . ", " .
            "      ".
            "$consSegment\n";
    }
    print "  " . $range->[4] . " " . sprintf("%7.2f",$range->[0]) . 
            " : " . sprintf("%$colWidths[0]d",$range->[1]) . " - " . sprintf("%$colWidths[1]d",$range->[2]) . 
            ", " . sprintf("%$colWidths[2]d",$msaToConPos[$range->[1]]) .  " - " . sprintf("%$colWidths[3]d",$msaToConPos[$range->[2]]) . 
            ", " . sprintf("%+4d",length($range->[3]) - ($msaToConPos[$range->[2]]-$msaToConPos[$range->[1]]+1)) . 
            ", " . " "x($msaToConPos[$range->[1]] - $msaToConPos[$minRange]) . $range->[3] . "\n";
  }
  print "\n";
  
}elsif ( $aggregationMethod eq "cluster" ) {
  # Cluster
  my $clusters = &clusterBlocks($mAlign, $ranges, 0, $copymin, $ratio, $cons, \%matrix);
  @tilingPath = @{$clusters};
  print "Clustered Ranges:\n";
  print "  Ratio   Range      Consensus\n";
  print "  -------------------------------\n";
  foreach my $range ( @tilingPath ) {
    print "" . sprintf("%7.2f",$range->[0]) . " : " . $range->[1] . " - " . $range->[2] . ", " . $range->[3] . "\n";
  }
  print "\n";
}

  # DEBUG:
  if ( 0 ) {
    my $blocks = [];
    foreach my $range ( @tilingPath ) {
      push @{$blocks}, [ $range->[1], $range->[2] ];
    }
    my $html = getMultAlnHTML(
      multAln => $mAlign,
      coloredBlocks => [
        ['#CD5C5C', $blocks],
      ],
      alignmentBlocks => $blocks,
      inclRef => 1,
      DEBUG => 1,
    );
    open OUT,">/home/rhubley/public_html/ri.html" or die;
    print OUT $html;
    close OUT;
  }


# Alter the consensus
#   - sort from high to low ranges keep range indexes valid during replacement
@tilingPath = sort { $a->[1] <=> $b->[1] } @tilingPath;
for ( my $i = $#tilingPath; $i >= 0; $i-- ) {
  my $range = $tilingPath[$i];
  my $oldCons = substr($cons, $range->[1], $range->[2]-$range->[1]+1);
  my $cOldCons = $oldCons;
  $cOldCons =~ s/-//g;

  print "Consensus Change:\n";
  print "" . printSequenceDiffs($cOldCons, $range->[3], "   ");
  print "\n";
 
  if ( $options{'interactive'} ) {
    print STDERR "s(kip),d(one) or press enter to keep\n";
    my $answer = <STDIN>;
    while ($answer !~ /^[sd]?$/ ) {
      print STDERR "Could not process $answer\nType 's','d' or press enter.\n";
      $answer = <STDIN>;
    }
    chomp $answer;
    last if ( $answer eq "d" );
    next if ( $answer eq "s" );
    substr($cons, $range->[1], $range->[2]-$range->[1]+1) = $range->[3];
  }else {
    substr($cons, $range->[1], $range->[2]-$range->[1]+1) = $range->[3];
  }
}
$cons =~ s/-//g;
print ">cons\n$cons\n";
print "\n";

if ( $options{'cons'} ) {
  open OUT,">$options{'cons'}" or die;
  print OUT ">cons\n$cons\n";
  close OUT;
}







################################################################################################

sub identifyUnstableBlockerWindows { 
  my $mAlign = shift;
  my $windowSizes = shift;
  my $copymin = shift;
  my $ratio = shift;
  my $matrix = shift;

  # create an array to translate consensus position to MSA position
  my $cons = $mAlign->consensus();
  my @consPos = ();
  my $consIdx = 0;
  for ( my $i = 0; $i < length($cons); $i++ ){
    my $cBase = substr($cons, $i, 1);
    next if ( $cBase eq "-" );
    $consPos[$consIdx] = $i;
    $consIdx++;
  }

  # variable length sliding window approach
  my @ranges = ();
  for ( my $i = 0; $i <= $#consPos; $i++ ) {
    foreach my $j ( @{$windowSizes} ) {
      if ( $i + $j < $#consPos ) {
        my $msaStart = $consPos[$i];
        my $msaEnd = $consPos[$i+$j-1];
        my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
                 &evalMSABlock( $mAlign, $msaStart, $msaEnd, $copymin, $ratio, $matrix, $cons);
        #print "BLOCKER: StartConsPos=$i, msa: $msaStart-$msaEnd WindowSize = $j Final: $numMajorityLen/$numOrigLen, $newseq, $newseqlen, $ncount\n";
        if ( $numMajorityLen != $numOrigLen ) {
          my $rratio;
          if ( $numOrigLen == 0 ) {
            $rratio = $numMajorityLen/0.1;
          }else {
            $rratio = $numMajorityLen/$numOrigLen;
          }
          push @ranges, [ $rratio, $msaStart, $msaEnd, $newseq ];   
          #print "StartConsPos=$i, WindowSize = $j Final: $numMajorityLen/$numOrigLen=$rratio , $newseq, $newseqlen, $ncount\n";
        }
      }
    }
  }
  return \@ranges;
}


# New clustering code
sub clusterBlocks {
  my $mAlign = shift;
  my $ranges = shift;
  my $gapAllowedDist = shift;
  my $copymin = shift;
  my $ratio = shift;
  my $cons = shift;
  my $matrix = shift;

  my @newRanges = ();
  # Sort by start column and longest block first
  my @currRanges = sort { $a->[1] <=> $b->[1] ||
                          $b->[2] <=> $a->[2] } @{$ranges};
  my $prevBlock;
  my $clusterStart;
  my $clusterEnd;
  my $blockIdx = 0;
  my @clusters = ();
  while ( $blockIdx <= $#currRanges ) {
    # Records = [ $rratio, $msaStart, $msaEnd, $newseq ];
    my $currBlock = $currRanges[$blockIdx];

    if ( ! defined $clusterStart ) {
      #print "    - Seeding cluster: " . $currBlock->[1] . " - " . $currBlock->[2] . " cluster_ratio = " . $currBlock->[0]  . "\n";
      $clusterStart = $currBlock->[1];
      $clusterEnd = $currBlock->[2];
      $blockIdx++;
      next;
    }

    # Check with overlap with current cluster
    if ( $currBlock->[1] - $gapAllowedDist - 1 <= $clusterEnd ) {
      # Evaluate the block to make sure it still has a high enough ratio
      my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
               &evalMSABlock( $mAlign, $clusterStart, $currBlock->[2], $copymin, $ratio, $matrix, $cons);
      if ( $numMajorityLen ) {
        #print "    - Appending to cluster: " . $currBlock->[1] . " - " . $currBlock->[2] . " cluster_ratio = " . ( $numMajorityLen/$numOrigLen ) . "\n";
        # Overlaps and has a significant ratio - expand the cluster
        $clusterEnd = $currBlock->[2];
        $blockIdx++;
        next;
      }else { 
        # Overlaps but has insignificant ratio
        # Emit cluster
        my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
              &evalMSABlock( $mAlign, $clusterStart, $clusterEnd, $copymin, $ratio, $matrix, $cons);
        push @clusters, [( $numMajorityLen/$numOrigLen ), $clusterStart, $clusterEnd, $newseq];
        #print "Cluster: $clusterStart - $clusterEnd " . ( $numMajorityLen/$numOrigLen ) . "\n";
        do { 
          $blockIdx++;
        }while ( $blockIdx <= $#currRanges && $currRanges[$blockIdx]->[1] <= $clusterEnd + $gapAllowedDist );
        undef $clusterStart;
        undef $clusterEnd;
        next;
      }
    }else {
      # Doesn't overlap cluster
      # Emit cluster
      my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
            &evalMSABlock( $mAlign, $clusterStart, $clusterEnd, $copymin, $ratio, $matrix, $cons);
      # TODO: With RuzzoTompa no guarantee this won't return a null result.
      #print "Eval cluster $clusterStart-$clusterEnd >$newseq< $numMajorityLen $numOrigLen\n";
      push @clusters, [( $numMajorityLen/$numOrigLen ), $clusterStart, $clusterEnd, $newseq];
      #print "Cluster: $clusterStart - $clusterEnd " . ( $numMajorityLen/$numOrigLen ) . "\n";
      do { 
        $blockIdx++;
      }while ( $blockIdx <= $#currRanges && $currRanges[$blockIdx]->[1] <= $clusterEnd + $gapAllowedDist );
      undef $clusterStart;
      undef $clusterEnd;
      next;
    }
  }
  if ( $clusterStart ) {
    my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
          &evalMSABlock( $mAlign, $clusterStart, $clusterEnd, $copymin, $ratio, $matrix, $cons);
    push @clusters, [( $numMajorityLen/$numOrigLen ), $clusterStart, $clusterEnd, $newseq];
    #print "Cluster: $clusterStart - $clusterEnd " . ( $numMajorityLen/$numOrigLen ) . "\n";
  }
  return \@clusters
}



#-----------------------------------------------------------------------------------------------------------------
##
##  Efficient Blocker evaluation
##
sub evalMSABlock{
  my $mAlign = shift;
  my $start = shift;
  my $end = shift;
  my $copymin = shift;
  my $ratio = shift;
  my $matrix = shift;
  my $consensus = shift;

  # This is an expensive operation.  If we are handed the
  # consensus, believe it to be derived from the mAlign object
  # previously and that the mAlign object hasn't been modified 
  # since that point.
  unless ( $consensus ) {
    $consensus = $mAlign->consensus();
  }

  # Get the alignment slice
  my ( $ref, $subMSA ) = $mAlign->getAlignmentBlock( start => $start, end => $end );
  #print "getAlignmentBlock( start=$start,end=$end) = $ref ". Dumper($subMSA) ."\n";
  # TODO: Do we always want to include the reference sequence?
  unshift @{$subMSA}, $ref;

  # Get the consensus slice
  $consensus = substr($consensus,$start,$end-$start+1);
  #print "Blocker:: consensus from $start-$end = $consensus\n";
  $consensus =~ s/-//g;
  my $conslen = length $consensus;

  # Calculate the sequence length distribution
  my $bestlength = 0;
  my $bestcnt = 0;
  my $secondbestlength = 0;
  my $secondbestcnt = 0;
  my $conscnt = 0;
  my $totalcnt = 0;
  my $maxlength = 0;
  my %cnt;
  foreach my $row ( @{$subMSA} ) {
    $totalcnt++;
    $row =~ s/-//g;
    my $seqlen = length($row);
    $conscnt++ if ( $seqlen == length($consensus) );
    $maxlength = $seqlen if ( $seqlen > $maxlength );
    $cnt{$seqlen}++;
    if ( $cnt{$seqlen} > $bestcnt ) {
      $bestlength = $seqlen;
      $bestcnt = $cnt{$seqlen};
    }elsif ( $cnt{$seqlen} > $secondbestcnt ) {
      $secondbestlength = $seqlen;
      $secondbestcnt = $cnt{$seqlen};
    }
  }
  #print "Best Length = $bestlength [ $bestcnt ]\n";
  #print "Second Best Length = $secondbestlength [ $secondbestcnt ]\n";
  #print "Count with consensus length = $conscnt\n";
  #print "Total rows " . ( $#{$subMSA}+1 ) . "\n";
  
  # Develop a crude consensus for each column
  my %col;
  foreach my $row ( @{$subMSA} ) {
    $row =~ s/-//g;
    if ( length($row) == $bestlength ) {
      my @string = split "", $row;
      for (my $i = 0; $i < $bestlength; ++$i) {
        $col{$i} .= $string[$i];
      }
    }
  }
  my $newconsensus = "";
  my $maxl;
  for (my $i = 0; $i < $bestlength; ++$i) {
    my @s = split "", $col{$i};
    my $max = -999999999;
    foreach my $let ( 'A','C','G','T','N') {
      my $score = 0;
      for (my $j = 0; $j <= $#s; ++$j) {
        $score += $matrix->{$let}{$s[$j]};
      }
      if ($score > $max) {
        $maxl = $let;
        $max = $score;
      }
    }
    $newconsensus .= "$maxl";
  }
  # TODO: Unclear why these were coded into AutoRunBlocker ( removed for now )
  #$newconsensus =~ s/^N+//; 
  #$newconsensus =~ s/N+$//;
  my $ncount = ( $newconsensus =~ tr/N/N/); 

  my $newlength = length($newconsensus);
  my $adjustedratio = $ratio;
  my $nratiocut = 10;
  if ($newlength > $conslen) {
    $adjustedratio -= 0.5;
    $nratiocut -= ($newlength - $conslen)/4;
  }
  if ( $bestcnt >= $copymin && (!$conscnt || $bestcnt/$conscnt >= $adjustedratio)
       && ($ncount < 2 || $newlength/$ncount > $nratiocut ) )
  {
    #print "ACCEPT: $start-$end bestcnt = $bestcnt ($copymin) conscnt = $conscnt adjustedratio = $adjustedratio ncount = $ncount  newlength = $newlength\n";
    return ( $bestcnt, $conscnt, $newconsensus, $newlength, $ncount, $secondbestcnt, $secondbestlength);
  }else {
    print "REJECT: $start-$end [" . ($#{$subMSA}+1) . " seqs in slice]  most_freq = $bestlength bp [$bestcnt] (must be >= $copymin), cons_len_freq = " . length($consensus) . " bp [$conscnt], most_freq/cons_freq = " . ($bestcnt/$conscnt) . " (must be >= $adjustedratio), ncount = $ncount (must be < 2) or newlength = $newlength (/ncount > $nratiocut),  new:$newconsensus\n";
    return ( 0, 0, "", 0, 0, 0, 0 );
  }
}



sub printSequenceDiffs {
  my $seq1 = uc(shift); # Old
  my $seq2 = uc(shift); # New
  my $prefix = shift;

  # For the purposes of labeling substitutions.
  my %mutChar = (
                  "CT" => 'i',
                  "TC" => 'i',
                  "AG" => 'i',
                  "GA" => 'i',
                  "GT" => 'v',
                  "TG" => 'v',
                  "GC" => 'v',
                  "CG" => 'v',
                  "CA" => 'v',
                  "AC" => 'v',
                  "AT" => 'v',
                  "TA" => 'v',
                  "NA" => "A",
                  "NC" => "C",
                  "NG" => "G",
                  "NT" => "T"
 );
  
  my $maxIdx = length($seq1);
  $maxIdx = length($seq2) if ( length($seq2) > $maxIdx );
  my $diffLine = "";
  for ( my $i = 0; $i < $maxIdx; $i++ ) {
    my $s1Base = " ";
    if ( $i < length($seq1) ) {
      $s1Base = substr($seq1, $i, 1);
    }
    my $s2Base = " ";
    if ( $i < length($seq2) ) {
      $s2Base = substr($seq2, $i, 1);
    }
    if ( exists $mutChar{$s1Base.$s2Base} ) {
      $diffLine .= $mutChar{$s1Base.$s2Base};
    }elsif ( $s2Base =~ /[BDHVRYKMSWNX]/ || 
             $s1Base =~ /[BDHVRYKMSWNX]/ ) {
      $diffLine .= "?";
    }elsif ( $s2Base eq " " || $s1Base eq " " ) {
      $diffLine .= "-";
    }else {
      $diffLine .= " ";
    }
  }
  return ( "    Orig: " . $seq1 . "\n" . "          ".$diffLine . "\n" . "     New: ".$seq2 . "\n" );
}
  
##
##  Original...Blocker 
##
sub Blocker{
  my $mAlign = shift;
  my $start = shift;
  my $end = shift;
  my $copymin = shift;
  my $ratio = shift;

  my $log;
  my ( $ref, $subMSA ) = $mAlign->getAlignmentBlock( start => $start, end => $end );
  #print "getAlignmentBlock( start=$start,end=$end) = $ref\n";
  my $consensus = $mAlign->consensus();
  $consensus = substr($consensus,$start,$end-$start+1);
  $log .= "con = $consensus\n";
  print "Blocker:: consensus from $start-$end = $consensus\n";
  $consensus =~ s/-//g;
  my $bestlength = 0;
  my $bestcnt = 0;
  my $secondbestlength = 0;
  my $secondbestcnt = 0;
  my $conscnt = 0;
  my $totalcnt = 0;
  my $maxsub = 0;
  my %cnt;
  unshift @{$subMSA}, $ref;
  foreach my $row ( @{$subMSA} ) {
    $log .= "$row\n";
    $totalcnt++;
    $row =~ s/-//g;
    my $seqlen = length($row);
    $conscnt++ if ( $seqlen == length($consensus) );
    $maxsub = $seqlen if ( $seqlen > $maxsub );
    $cnt{$seqlen}++;
    if ( $cnt{$seqlen} > $bestcnt ) {
      $bestlength = $seqlen;
      $bestcnt = $cnt{$seqlen};
    }elsif ( $cnt{$seqlen} > $secondbestcnt ) {
      $secondbestlength = $seqlen;
      $secondbestcnt = $cnt{$seqlen};
    }
  }
  $log .= "Best Length = $bestlength [ $bestcnt ]\n";
  $log .= "Second Best Length = $secondbestlength [ $secondbestcnt ]\n";
print "Best Length = $bestlength [ $bestcnt ]\n";
print "Second Best Length = $secondbestlength [ $secondbestcnt ]\n";
  
  # For the purposes of labeling substitutions.
  my %mutChar = (
                  "CT" => 'i',
                  "TC" => 'i',
                  "AG" => 'i',
                  "GA" => 'i',
                  "GT" => 'v',
                  "TG" => 'v',
                  "GC" => 'v',
                  "CG" => 'v',
                  "CA" => 'v',
                  "AC" => 'v',
                  "AT" => 'v',
                  "TA" => 'v' );
  
  my $conslen = length $consensus;
  
  # to report a cluster of longer sequences; these tend to be off the screen
  my $outlier = "";
  if ($maxsub > $bestlength + 3) {
    my $diffinlen = 0;
    my $diffinlen2 = $secondbestlength - $bestlength;
    for (my $i = $bestlength + 3; $i <= $maxsub; $i++) {
      my $lengthdiff = $i - $bestlength;
      if ($cnt{$i} && $cnt{$i} > 3 && $cnt{$i} > $bestcnt/10 &&
        $i - $secondbestlength > 2 && $lengthdiff > $diffinlen2) {
        $diffinlen = $i - $bestlength;
        $outlier = ", $cnt{$i} copies of $i bp";
      }
    }
  }
  
  my %col;
  foreach my $row ( @{$subMSA} ) {
    $row =~ s/-//g;
    if ( length($row) == $bestlength ) {
      my @string = split "", $row;
      for (my $i = 0; $i < $bestlength; ++$i) {
        $col{$i} .= $string[$i];
      }
    }
  }

  open (IN, "$FindBin::RealBin/../Matrices/linupmatrix ") or die;
  my @cons;
  my @l;
  my %matrix;
  while (<IN>) {
    s/^\s+//;
    if (/^A/) {
      @l = split;
      @cons = @l;
    } elsif (/^[\d-]/) {
      my @s = split;
      my $l = shift @cons;
      for (my $i = 0; $i <= $#s; ++$i) {
        $matrix{$l}{$l[$i]} = $s[$i];
      }
    }
  }
  close IN;
  
  my $printalert = "";
  $printalert = "Length difference. " unless $bestlength == (length $consensus);
  my @oldstring = split "", $consensus;
  my $newstring = "";
  my $maxl;
  my $spacer = "    ";
  for (my $i = 0; $i < $bestlength; ++$i) {
    my @s = split "", $col{$i};
    my $max = -999999999;
    foreach my $let ( 'A','C','G','T','N') {
      my $score = 0;
      for (my $j = 0; $j <= $#s; ++$j) {
        $score += $matrix{$let}{$s[$j]};
      }
      if ($score > $max) {
        $maxl = $let;
        $max = $score;
      }
    }
    $newstring .= "$maxl";
    if ($i <= $#oldstring) {
      if ($maxl ne 'N' && $oldstring[$i] eq 'N' ||
  	$maxl =~ /[AG]/ && $oldstring[$i] =~ /[CT]/ ||
  	$maxl =~ /[CT]/ && $oldstring[$i] =~ /[AG]/) {
        $printalert = "Sequence difference. " unless $printalert;
        $spacer .= $maxl;
      } else {
        $spacer .= " ";
      }
    }
  }
  
  ## Calculate a visual diff line between the old/new consensus
  my $diff = $spacer;
  if ( $options{'cmdiffs'} ) {
    $diff = "    ";
    my $diffLen = length($consensus);
    $diffLen = length($newstring) if ( length($newstring) > $diffLen );
    for ( my $i = 0; $i < $diffLen; $i++ )
    {
      my $newB;
      $newB = substr( $newstring, $i, 1 ) if ( $i < length($newstring) );
      my $oldB;
      $oldB = substr( $consensus, $i, 1 ) if ( $i < length($consensus));
      if ( ! defined $newB || ! defined $oldB ) {
        $diff .= "-";
      }elsif ( my $mc = $mutChar{ uc( $newB . $oldB ) } )
      {
        $diff .= $mc;
      }elsif (    ( $newB =~ /[BDHVRYKMSWNX]/i )
               || ( $oldB =~ /[BDHVRYKMSWNX]/i ) )
      {
        $diff .= "?";
      }else
      {
        $diff .= " ";
      }
    }
  }
  
  if ($printalert) {
    if ($printalert =~ /^Len/) {
      $log .= "\n$printalert$bestcnt of $totalcnt for $bestlength bp ($secondbestcnt second best for $secondbestlength bp$outlier, $conscnt for original $conslen bp)\nold $consensus\n$diff\nnew $newstring\n";
    } else {
      $log .= "\n$printalert$bestcnt of $totalcnt for $bestlength bp ($secondbestcnt second best for $secondbestlength bp$outlier)\nold $consensus\n$diff\nnew $newstring\n";
    }
  } else {
    $log .= "\n$bestcnt of $totalcnt for $bestlength bp ($secondbestcnt second best for $secondbestlength bp$outlier)\nsame $consensus\n";
  }

  $newstring =~ s/^N+//; 
  $newstring =~ s/N+$//;
  my $ncount = ( $newstring =~ tr/N/N/); 

  my $newlength = length($newstring);
  my $adjustedratio = $ratio;
  my $nratiocut = 10;
  if ($newlength > length($consensus)) {
    $adjustedratio -= 0.5;
    $nratiocut -= ($newlength - length($consensus))/4;
  }
  if ( $bestcnt >= $copymin && (!$conscnt || $bestcnt/$conscnt >= $adjustedratio)
       && ($ncount < 2 || $newlength/$ncount > $nratiocut ) )
  {
    return ( $bestcnt, $conscnt, $newstring, $newlength, $ncount, $secondbestcnt, $secondbestlength, $log);
  }else {
    print "bestcnt = $bestcnt ($copymin) conscnt = $conscnt ncount = $ncount  newlength = $newlength\n";
    return ( 0, 0, "", 0, 0, 0, 0 );
  }
}


## New code for visualizing MultAln object in html
## TODO: return the offset to the first position of the MSA for aligning more information
## 
sub getMultAlnHTML {
    my %parameters = @_;

    my $mAlign  = $parameters{'multAln'}
      or croak "getMultAlnHTML() missing 'multAln' parameter!\n";

    # If not provided, assume 'inclRef' is false
    my $inclRef = $parameters{'inclRef'} ? 1 : 0;
    my $DEBUG   = $parameters{'DEBUG'}   || 0;

    # We'll store all HTML in this scalar and return it
    my $html = '';

    #-----------------------------------------------------------------------------
    # 1) Generate HTML header / doctype
    #-----------------------------------------------------------------------------
    $html .= <<"END_HEADER";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
        "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <title>Alignments</title>
  <style type="text/css">
    /* Example CSS classes -- not strictly needed if we do inline styles for color blocks */
    font.lowQual {
      background: #CD5C5C;
    }
    font.deletion {
      background: #FF0000;
    }
    font.duplication {
      background: #0000FF;
    }
    font.unknown {
      background: #FFFF00;
    }
    font.dupFlank {
      background: #C0C0C0;
    }
  </style>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
</head>
<body>
<PRE>
END_HEADER

    #-----------------------------------------------------------------------------
    # 2) Compute left/right flank padding for alignment lines
    #-----------------------------------------------------------------------------
    my ($maxLeftLen, $maxRightLen) = _compute_flank_padding($mAlign);

    #-----------------------------------------------------------------------------
    # 3) Print optional MSA column header (the top numeric scale)
    #    - Often you do this with _build_column_header, or skip if you want none.
    #-----------------------------------------------------------------------------
    my $alnLength = length($mAlign->getReferenceSeq());
    $html .= _build_column_header($alnLength, $maxLeftLen);

    #-----------------------------------------------------------------------------
    # 4) Print any custom "vertical headers" the user asked for
    #    e.g. [ ["lowQualScore",[0,0,1.2,...]], ["GCcontent",[0.4,0.7,0.0,...]] ]
    #-----------------------------------------------------------------------------
    if (my $headers = $parameters{'headersAboveMSA'}) {
        $html .= _print_vertical_headers($headers, $maxLeftLen);
    }

    #-----------------------------------------------------------------------------
    # 5) Print the CONSENSUS line, optionally including reference in consensus
    #-----------------------------------------------------------------------------
    {
        my $consensus = $inclRef 
            ? $mAlign->consensus(inclRef => 1)
            : $mAlign->consensus();

        my $name       = "consensus";
        my $namePad    = 30 - length($name);
        my $label      = "<b><i>$name</i></b>" . (' ' x $namePad) . ": ";
        my $leftSpaces = ' ' x $maxLeftLen;

        # Annotate each base with a hover indicating ungapped position
        my $annotatedSeq = _annotate_positions_for_hover($consensus, "Cons pos");

        $html .= $label 
              .  $leftSpaces
              .  qq(<font color="blue">$annotatedSeq</font>\n);
    }

    #-----------------------------------------------------------------------------
    # 6) Print the REFERENCE line, also annotated for hover if desired
    #-----------------------------------------------------------------------------
    {
        my $refName    = "Reference ( " . $mAlign->getReferenceName() . " )";
        my $pad        = 30 - length($refName);
        my $label      = "<b><i>$refName</i></b>" . (' ' x $pad) . ": ";
        my $leftSpaces = ' ' x $maxLeftLen;

        my $refSeq     = $mAlign->getReferenceSeq();
        my $annotated  = _annotate_positions_for_hover($refSeq, "Ref pos");

        $html .= $label
              .  $leftSpaces
              .  qq(<font color="blue">$annotated</font>\n);
    }

    #-----------------------------------------------------------------------------
    # 7) Build a color array for each aligned sequence if 'coloredBlocks' is given
    #    Format:  coloredBlocks => [
    #       [ '#FF0000', [ [start,end], [start,end], ... ] ],
    #       [ '#00FF00', [ [start,end] ] ],
    #       ...
    #    ]
    #-----------------------------------------------------------------------------
    my $coloredBlocksRef = $parameters{'coloredBlocks'} || [];
    # We'll combine them into a single array of "per-column color" for each sequence on the fly.

    #-----------------------------------------------------------------------------
    # 8) Print each ALIGNED sequence
    #-----------------------------------------------------------------------------
    for (my $i = 0; $i < $mAlign->getNumAlignedSeqs(); $i++) {

        my $seqName = $mAlign->getAlignedName($i);
        $seqName    = substr($seqName, 0, 30) if length($seqName) > 30;
        my $pad     = 30 - length($seqName);
        my $label   = "<b>$seqName</b>" . (' ' x $pad) . ": ";

        my $lfSeq   = $mAlign->getLeftFlankingSequence($i);
        my $alnSeq  = $mAlign->getAlignedSeq($i);
        my $rfSeq   = $mAlign->getRightFlankingSequence($i);

        # Spaces for left flank
        my $lfShift = length($lfSeq) - $mAlign->getAlignedStart($i);
        my $padding = ' ' x ($maxLeftLen - $lfShift);

        # Possibly highlight left flank in blue
        my $leftFlank = (defined $parameters{'leftFlankingID'}
                        && $parameters{'leftFlankingID'} eq ($i - 1))
                        ? qq(<font color="blue">) . lc($lfSeq) . "</font>"
                        : lc($lfSeq);

        # Possibly highlight right flank in blue
        my $rightFlank = (defined $parameters{'rightFlankingID'}
                         && $parameters{'rightFlankingID'} eq ($i - 1))
                         ? qq(<font color="blue">) . lc($rfSeq) . "</font>"
                         : lc($rfSeq);

        # Apply color blocks to the interior alignment
        my $coloredSeq = _apply_colored_blocks(
            $alnSeq, 
            $mAlign->getAlignedStart($i),
            $coloredBlocksRef
        );

        # Combine them all into one line
        my $line = $label
                 . $padding
                 . $leftFlank
                 . $coloredSeq
                 . $rightFlank
                 . "\n";
        $html .= $line;
    }

    #-----------------------------------------------------------------------------
    # 9) If we have alignmentBlocks => [ [start,end], [start,end], ... ],
    #    print them at the bottom. (like your old _print_histogram_info)
    #-----------------------------------------------------------------------------
    if (my $blocks = $parameters{'alignmentBlocks'}) {
        # Only call if we actually have an array of [start,end]
        if (@$blocks) {
            $html .= "\n\n";
            $html .= _print_alignment_blocks($mAlign, $blocks, $maxLeftLen);
        }
    }

    #-----------------------------------------------------------------------------
    # Finish and return
    #-----------------------------------------------------------------------------
    $html .= "</PRE>\n</body>\n</html>\n";
    return $html;
}

###############################################################################
# Helper Subroutines
###############################################################################

#------------------------------------------------------------------------------
# _compute_flank_padding
#------------------------------------------------------------------------------
sub _compute_flank_padding {
    my ($mAlign) = @_;
    my $maxLeftLen  = 0;
    my $maxRightLen = 0;

    my $maxQueryEnd = length($mAlign->getReferenceSeq());
    for my $seqNum (0 .. $mAlign->getNumAlignedSeqs() - 1) {
        my $relLeftLen = length($mAlign->getLeftFlankingSequence($seqNum))
                       - $mAlign->getAlignedStart($seqNum);
        my $relRightLen = length($mAlign->getRightFlankingSequence($seqNum))
                        - ($maxQueryEnd - $mAlign->getAlignedEnd($seqNum));
        $maxLeftLen  = $relLeftLen  if $relLeftLen  > $maxLeftLen;
        $maxRightLen = $relRightLen if $relRightLen > $maxRightLen;
    }
    return ($maxLeftLen, $maxRightLen);
}

#------------------------------------------------------------------------------
# _build_column_header
#------------------------------------------------------------------------------
# If you still want a column header for 1..N columns, do so vertically here.
sub _build_column_header {
    my ($alnLength, $maxLeftLen) = @_;
    return '' if $alnLength < 1;

    my $html      = '';
    my $maxDigits = length($alnLength);

    # Pre-build an array of empty strings, one per digit-row
    my @rows = map { '' } (1..$maxDigits);

    # For each column number, create a right-aligned string of length $maxDigits
    for my $col (1 .. $alnLength) {
        my $colStr = sprintf("%*d", $maxDigits, $col);  # e.g. '  42' if maxDigits=3
        # Distribute each digit into the corresponding row
        for my $d (0 .. $maxDigits - 1) {
            $rows[$d] .= substr($colStr, $d, 1);
        }
    }

    # Print each row, preceded by 30 spaces for the "label" and $maxLeftLen
    for my $r (0 .. $maxDigits - 1) {
        my $namePad = "MSA Column" . ' ' x 20;
        my $lfPad   = ' ' x $maxLeftLen;
        $html .= $namePad . ": " . $lfPad . $rows[$r] . "\n";
    }
    $html .= "\n";  # Blank line
    return $html;
}

#------------------------------------------------------------------------------
# _annotate_positions_for_hover
#------------------------------------------------------------------------------
# For each non-gap character, wrap in <span title="..."> so you can see
# the ungapped position on mouseover.
sub _annotate_positions_for_hover {
    my ($sequence, $labelPrefix) = @_;
    my $html         = '';
    my $ungappedPos  = 1;

    for (my $i = 0; $i < length($sequence); $i++) {
        my $char = substr($sequence, $i, 1);
        if ($char eq '-') {
            # Gap: do not increment pos
            $html .= $char;
        }
        else {
            my $title = qq(title="$labelPrefix: $ungappedPos");
            $html .= qq(<span $title>$char</span>);
            $ungappedPos++;
        }
    }
    return $html;
}

#------------------------------------------------------------------------------
# _print_vertical_headers
#------------------------------------------------------------------------------
# Given an array of:
#   [ [ $title, [ colVal1, colVal2, ... colValN ] ],
#     [ $title2, [ colVal1, colVal2, ... colValN ] ], ... ]
# For each "header," we print a vertical stacking of digits (like your
# prior "score array" code). Zero values become blank.
sub _print_vertical_headers {
    my ($headers, $maxLeftLen) = @_;
    my $html = '';

    foreach my $header (@$headers) {
        my ($title, $colVals) = @$header;
        next if !$colVals || !@$colVals;  # skip empty arrays

        # 1) Convert each val to a string
        my $maxWidth = 0;
        my @formatted;
        for my $val (@$colVals) {
            # If exactly 0, we make it blank
            my $str = ($val == 0) ? '' : sprintf("%0.1f", $val);
            $maxWidth = length($str) if length($str) > $maxWidth;
            push @formatted, $str;
        }

        # 2) We'll have $maxWidth rows
        my @lines = map { '' } (1..$maxWidth);

        for my $str (@formatted) {
            my $padded = (' ' x ($maxWidth - length($str))) . $str;
            for (my $i = 0; $i < $maxWidth; $i++) {
                $lines[$i] .= substr($padded, $i, 1);
            }
        }

        # 3) Print them
        my $paddedLabel = $title . ' ' x (30 - length($title));
        for my $line (@lines) {
            $html .= $paddedLabel . ": " . (' ' x $maxLeftLen) . $line . "\n";
        }
        $html .= "\n";
    }

    return $html;
}

#------------------------------------------------------------------------------
# _apply_colored_blocks
#------------------------------------------------------------------------------
# We build an array @colorForPos that tells us what color each column
# in $alignedSeq should have (if any). Then we walk the sequence
# to build an HTML string with <span style="background-color:...">.
sub _apply_colored_blocks {
    my ($alignedSeq, $alignedStart, $coloredBlocksRef) = @_;

    my $seqLen = length($alignedSeq);
    # Initialize no color for each position
    my @colorForPos = (undef) x $seqLen;

    # For each [ color, [ [start,end], [start,end], ... ] ] tuple
    # fill colorForPos for the overlapping columns
    foreach my $blockDef (@$coloredBlocksRef) {
        my ($color, $ranges) = @$blockDef;
        next if !$ranges || !@$ranges;

        # For each [start,end]
        foreach my $range (@$ranges) {
            my ($start, $end) = @$range;
            # Convert them to local coords for this aligned sequence
            my $localStart = $start - $alignedStart;
            my $localEnd   = $end   - $alignedStart;

            # Skip if out of range
            next if $localEnd < 0 || $localStart >= $seqLen;

            # Clamp
            $localStart = 0          if $localStart < 0;
            $localEnd   = $seqLen-1  if $localEnd   >= $seqLen;

            # Fill in color
            for my $pos ($localStart .. $localEnd) {
                $colorForPos[$pos] = $color;
            }
        }
    }

    # Now build the final HTML string
    my $html  = "<b>";
    my $pos   = 0;

    while ($pos < $seqLen) {
        my $currentColor = $colorForPos[$pos];
        my $char         = substr($alignedSeq, $pos, 1);

        # If there's a color here, we start a <span>; otherwise no color
        if (defined $currentColor) {
            # Find how many consecutive columns share this color
            my $runStart = $pos;
            while ($pos < $seqLen && defined $colorForPos[$pos] && $colorForPos[$pos] eq $currentColor) {
                $pos++;
            }
            # $pos is now the first position after this color run
            my $runStr = substr($alignedSeq, $runStart, $pos - $runStart);
            $html .= qq(<span style="background-color:$currentColor">$runStr</span>);
        }
        else {
            # No color here
            $html .= $char;
            $pos++;
        }
    }

    $html .= "</b>";
    return $html;
}

#------------------------------------------------------------------------------
# _print_alignment_blocks
#------------------------------------------------------------------------------
# If the user wants to list each [start,end] block at the bottom, e.g., to show
# the alignment sequences that appear in that region, or any summary lines, etc.
# For demonstration, we'll do something simple. Adapt as needed.
sub _print_alignment_blocks {
    my ($mAlign, $blocks, $maxLeftLen) = @_;
    return '' if !$blocks || !@$blocks;

    # Example: collect sequences from each block and display
    # (Similar to older "histogram" approach)
    my $html = '';
    my $label        = "blockSeqs";
    my $paddedLabel  = $label . " " x (30 - length($label));

    my $maxRows = 0;
    my @columnSeqs;
    foreach my $block (@$blocks) {
        my ($start, $end) = @$block;
        my ($cons, $seqsRef) = $mAlign->getAlignmentBlock(
            start        => $start,
            end          => $end,
            rawSequences => 1
        );

        # Sort them by length descending
        my @sorted = sort { length($b) <=> length($a) } @$seqsRef;
        push @columnSeqs, \@sorted;                      
        $maxRows = @sorted if @sorted > $maxRows;           
    }

    my $label        = "blockSeqs";                                               
    my $paddedLabel  = $label . " " x (30 - length($label));                
    # Build each row                                                                               
    for my $rowIdx (0 .. $maxRows - 1) {                             
        my $line = ' ' x $maxLeftLen;                 
        my $pos  = 0;

        for my $colIdx (0 .. $#$blocks) {
            my ($start, $end) = @{$blocks->[$colIdx]};
            my $colWidth = $end - $start + 1;

            # Insert spaces for any gap before this block
            if (($start - $pos) > 0) {
                $line .= ' ' x ($start - $pos);
            }

            # Grab next seq from column
            my $colSeqArray = $columnSeqs[$colIdx];
            my $colSeq      = ($rowIdx < @$colSeqArray) ? $colSeqArray->[$rowIdx] : '';
            $colSeq         = '.' if $colSeq eq '';

            # Pad to fill the block width
            if (length($colSeq) < $colWidth) {
                $colSeq .= ' ' x ($colWidth - length($colSeq));
            }
            $line .= $colSeq;
            $pos   = $end + 1;
        }
        $html .= "$paddedLabel: $line\n";
    }

    $html .= "\n\n";
    return $html;
}

