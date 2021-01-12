#!/usr/bin/perl -w
###---------------------------------------------------------------------------##
###  File:
###      @(#) extendFlankingSeqs.pl
###  Author:
###      Robert M. Hubley   rhubley@systemsbiology.org
###  Description:
###      Read in a crossmatch out file and extend the sequences.
###
##******************************************************************************
##*  This software is provided ``AS IS'' and any express or implied            *
##*  warranties, including, but not limited to, the implied warranties of      *
##*  merchantability and fitness for a particular purpose, are disclaimed.     *
##*  In no event shall the authors or the Institute for Systems Biology        *
##*  liable for any direct, indirect, incidental, special, exemplary, or       *
##*  consequential damages (including, but not limited to, procurement of      *
##*  substitute goods or services; loss of use, data, or profits; or           *
##*  business interruption) however caused and on any theory of liability,     *
##*  whether in contract, strict liability, or tort (including negligence      *
##*  or otherwise) arising in any way out of the use of this software, even    *
##*  if advised of the possibility of such damage.                             *
##*                                                                            *
##******************************************************************************

=head1 NAME

 extendFlankingSeqs.pl - Extend sequences ranges defined by crossmatch out file

=head1 SYNOPSIS

 extendFlankingSeqs.pl [-version] [-h(elp)] 
                       [[-r(ight_flank) #] [-l(eft_flank) #] | [-f(lank) #]]
                       [-g(ap) #] [-re(index)] [-s(hrink)]
                       -d(atabase) <2bit file> 
                       -i(nput) <cross_match file> -o(utput) <fasta file>

=head1 DESCRIPTION

Given an set of hits ( *.out or *.align ) in crossmatch format of a set of genomic
sequences against a reference (consensus) sequence, add flanking sequence to the aligned
ranges and output in FASTA format.  The alignment is expected to be query=genomic and
subject=reference. E.g:

  cross_match genomic_copies.fa consensus.fa 
or
  rmblast.pl genomic_copies.fa consensus.fa

The amount of base pairs to add to the aligned sequences may be specified for both sides
using the -flank parameter or seperately as -left_flank/-right_flank.  The left/right
sides refer to the 5'/3' of the reference sequence.  

The sequences are pulled from a 2bit database file and a index is created of that file
the first time the script is on the database.  

The options are:

=over 4

=item -version

Displays the version of the program

=item -right_flank #

The number of bp to extend the right side of the aligned sequences.  Right/Left
refers to the side of the reference sequence used in the alignment.  This script 
assumes that the subject sequence was the reference and the genomic sequences were 
used as the query.

=item -left_flank #

The number of bp to extend the left side of the aligned sequences.  Right/Left
refers to the side of the reference sequence used in the alignment.  This script 
assumes that the subject sequence was the reference and the genomic sequences were 
used as the query.

=item -flank #

The number of bp to extend both right/left sides of the aligned sequences.

=item -gap #

Overlapping ranges are automatically merged, however it is sometimes advantageous
to merge nearby ranges as well.  If a -gap value is specified, ranges separated by
the specified amount will also be merged.

=item -shrink

In some cases the alignments do not fully extend to the edges of the source 
sequences.  In this case the program doesn't need to extend the sequences and
by default leaves the original ranges intact.  If -shrink is specified the
ranges can also be shrunk to be the alignment ranges + the specified flanking
number of bases.

=item -reindex

An index file is created from the 2bit database if one doesn't exist.  Each time
the program is run on the same database the same index file will be used to speed
up repeated runs.  If this option is specified it will force rebuilding of this
index file.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2019-2021 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Storable qw(nstore retrieve);
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use NCBIBlastSearchEngine;
use SearchResult;
use SearchResultCollection;

# Program version
my $Version = $RepModelConfig::VERSION;

#
# Paths
#
my $config           = $RepModelConfig::configuration;
my $TWOBITINFO       = $config->{'UCSCTOOLS_DIR'}->{'value'} . "/twoBitInfo";
my $TWOBITTOFA       = $config->{'UCSCTOOLS_DIR'}->{'value'} . "/twoBitToFa";
my $FATOTWOBIT       = $config->{'UCSCTOOLS_DIR'}->{'value'} . "/faToTowBit";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-version',
                    '-reindex|re',
                    '-shrink|s',
                    '-gap|g=i',
                    '-input|i=s',
                    '-output|o=s',
                    '-database|d=s',
                    '-left_flank|l=i',
                    '-right_flank|r=i',
                    '-flank|f=i',
);
my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0 | less";
  exit( 1 );
}

if ( ! exists $options{'database'} || ! -s $options{'database'} ) {
  die "\n\nDatabase doesn't exist!  Must supply a TwoBit file!\nUse '$0 -h' to view the help.\n\n";
}

if ( ! exists $options{'input'} || ! -s $options{'input'} ){
  die "\n\nMust supply input alignments in crossmatch format!\nUse '$0 -h' to view the help.\n\n";
}

if ( ! exists $options{'output'} ){
  die "\n\nMust supply output filename!\nUse '$0 -h' to view the help.\n\n";
}

my $inputFile = $options{'input'};
my $outputFile = $options{'output'};

my $leftFlank = 100;
my $rightFlank = 100;

if ( $options{'flank'} ) {
  $leftFlank = $rightFlank = $options{'flank'};
}else {
  $leftFlank = $options{'left_flank'} if ( $options{'left_flank'} );
  $rightFlank = $options{'right_flank'} if ( $options{'right_flank'} );
}

my $gapTolerance = 0;
$gapTolerance = $options{'gap'} if ( $options{'gap'} );

my $twoBitFile = $options{'database'};
my($twoBitFilename, $twoBitDir, $twoBitSuffix) = fileparse($twoBitFile);
my $indexFile = "$twoBitFilename.idx";

if ( ! -e $indexFile || exists $options{'reindex'} ) {
  # Must rebuild master index
  my %seqLens = ();
  open IN,"$TWOBITINFO $twoBitFile stdout|" or die "Could not run $TWOBITINFO on $twoBitFile!\n";
  my $lines = 0;
  while ( <IN> ) {
    $lines++;
    if ( /^(\S+)\s+(\d+)/ ) {
      $seqLens{$1} = $2;
    }
  }
  close IN;
  nstore(\%seqLens,$indexFile);
}

if ( ! -s $indexFile ) {
  die "Index file missing or empty!\n";
}
my $seqLens = retrieve($indexFile);

&addFlankingSequence($twoBitFile, $inputFile, $outputFile, $seqLens, $leftFlank, $rightFlank);


#--------------------------------------------------------------------------------------

## 
##  Pulled in addFlankingSequence.pl to improve performance.  The
##  original program needs to read in the details of the *.2bit file
##  each time to get the sequence lengths.  By pulling this into
##  the program which calls it many times, we can perform this 
##  function in advance and keep the sequence size datastructure in
##  memory.
##
##  Given a *.align or *.out file containing alignments of individual
##  genomic sequences to a consensus, grab sequences the squences
##  the genome and ensure each has at least X bp flanking unaligned
##  sequence attached.  
##
##  The assumption is that the individual query sequences are already
##  isolated from the genome as subsequences and may contain 0 or
##  more flanking bases.  The nomenclature for these query sequences
##  MUST be id_start_end or id_start_end_R ( for reverse orientation ).
##  For now it is assumed that these start/end positions are 1-based
##  fully-closed coordinates.
##
sub addFlankingSequence {
  my $twoBitFile = shift;
  my $alignFile = shift;
  my $outputFile = shift;
  my $twoBitSeqSizes = shift;
  my $flankLeftSize = shift;
  my $flankRightSize = shift;
  
  # Renaming to be compatible with original external script
  # TODO: clean up
  my $genome    = $twoBitFile;
  my $flankleft = 100;
  $flankleft = $flankLeftSize if ( defined $flankLeftSize );
  my $flankright = 100;
  $flankright = $flankRightSize if ( defined $flankRightSize );

  if ( ! -s $genome ) {
    die "addFlankingSequence.pl: Could not find genome file $genome\n";
  }

  my $DEBUG = 0;
  my %seqLens = %{$twoBitSeqSizes};
 
  my %ranges    = ();
  my $leftAdequate = 0;
  my $rightAdequate = 0;
  my $leftUnextendable = 0;
  my $rightUnextendable = 0;
  my $totalRanges = 0;
  open IN,"<$alignFile" or die "Could not open $alignFile for reading!\n";
  while ( <IN> )
  {
    # Standard cross_match-like aligments of repseq (query) vs rep (consensus) as is created
    # by dothemsimple.pl.
    #   2514  8.77 1.46 0.58  DS497995.1:134929-141580_R     6360  6701 (50)  C rnd-2_family-1#LINE/L1   (10)  6866  6522
    #   33545 6.47 0.57 0.02  CM000382.2:27738956-27743162_R  53  4256 (50)    rnd-2_family-1#LINE/L1   2132  6358 (518)  
    if ( /^\s*(\d+\s+\d+\.\d+.*)$/ ) 
    {
      my @flds = split;
      my $raw_id = $flds[4];
      my $align_start  = $flds[5];
      my $align_end    = $flds[6];
      my $align_remain = $flds[7];
      $align_remain =~ s/[\)\(]//g;
      my $align_orient = "+";
      $align_orient = "-" if ( $flds[8] eq "C" );
  
      my $seqRecordID = "";
      my $seqRecordStart;
      my $seqRecordEnd;
      my $seqRecordOrient;
      my $seqRecordLen;
      # These have traditionally been in 1-based full-closed coordinates
      if ( $raw_id =~ /^(\S+)\_(\d+)\_(\d+)\_?(R?)$/ ||
           $raw_id =~ /^(\S+):(\d+)-(\d+)/ )
      {
        $seqRecordID = $1;
        $seqRecordStart = $2;
        $seqRecordEnd = $3;
        $seqRecordOrient = "+";
        # Keep ranges in low to high order
        if ( $seqRecordStart > $seqRecordEnd ) {
          $seqRecordStart = $3;
          $seqRecordEnd = $2;
          $seqRecordOrient = "-";
        }elsif ( defined $4 &&  $4 eq "R" ) {
          $seqRecordOrient = "-";
        }
        die "Could not find $seqRecordID in sequence index!\n" if ( !exists $seqLens{$seqRecordID} );
        $seqRecordLen = $seqLens{$seqRecordID};
      }else {
        $seqRecordID = $raw_id;
        $seqRecordOrient = "+";
        $seqRecordStart = 1;
        die "Could not find $seqRecordID in sequence index!\n" if ( !exists $seqLens{$seqRecordID} );
        $seqRecordEnd = $seqRecordLen = $seqLens{$seqRecordID};
      }

      #print "SeqRecord: $seqRecordID:$seqRecordStart-$seqRecordEnd ($seqRecordOrient) <=> Alignment: $align_start-$align_end [$align_remain] ($align_orient)\n";

      my $globalStart;
      my $globalEnd;

      # SeqRecord Alignment
      # --------- ---------
      #     +        +          globalStart=left=seqRecordStart+AlignmentStart-1, globalEnd=right=seqRecordStart+AlignmentEnd-1 [+]
      #     +        -          globalStart=right=seqRecordStart+AlignmentStart-1, globalEnd=left=seqRecordStart+AlignmentEnd-1 [-]
      #     -        +          globalStart=right=seqRecordEnd-AlignmentEnd+1, globalEnd=left=seqRecordEnd-AlignmentStart+1     [+]
      #     -        -          globalStart=left=seqRecordEnd-AlignmentEnd+1, globalEnd=right=seqRecordEnd-AlignmentStart+1     [-]
      if ( $seqRecordOrient eq "+" ) {
        # forward seqRecord
        $globalStart = $seqRecordStart + $align_start - 1;
        $globalEnd = $seqRecordStart + $align_end - 1;
      }else {   
        # reversed seqRecord 
        $globalStart = $seqRecordEnd - $align_end + 1;
        $globalEnd = $seqRecordEnd - $align_start + 1;
      }
      #print "Adjusted1: $globalStart-$globalEnd\n";
        
      if ( $align_orient eq $seqRecordOrient ) {
        $leftUnextendable++ if ( $globalStart == 1 );  
        $rightUnextendable++ if ( $globalEnd == $seqRecordLen );  
        $globalStart -= $flankleft;
        $globalStart = 1 if ( $globalStart < 1 );
        $globalEnd += $flankright;
        $globalEnd = $seqRecordLen if ( $globalEnd > $seqRecordLen );
      }else {
        $rightUnextendable++ if ( $globalStart == 1 );  
        $leftUnextendable++ if ( $globalEnd == $seqRecordLen );  
        $globalStart -= $flankright;
        $globalStart = 1 if ( $globalStart < 1 );
        $globalEnd += $flankleft;
        $globalEnd = $seqRecordLen if ( $globalEnd > $seqRecordLen );
      }
      #print "Adjusted2: $globalStart-$globalEnd\n";
 
      if ( $align_orient eq "+" ) {
        if ( $align_start-1 >= $flankleft ){
          $leftAdequate++;  
          unless( $options{'shrink'} ) {
            $globalStart = $seqRecordStart;
          }
        }
        if ( $align_remain >= $flankright ){
          $rightAdequate++;  
          unless ( $options{'shrink'} ) {
            $globalEnd = $seqRecordEnd;
          }
        }
      }else {
        if ( $align_start-1 >= $flankright ){
          $rightAdequate++;  
          unless ( $options{'shrink'} ) {
            $globalEnd = $seqRecordEnd;
          }
        }
        if ( $align_remain >= $flankleft ){
          $leftAdequate++;  
          unless ( $options{'shrink'} ) {
            $globalStart = $seqRecordStart;
          }
        }
      }
     
      #print " Range: $globalStart-$globalEnd ($align_orient)\n";
      push @{$ranges{$seqRecordID}}, [$globalStart, $globalEnd, $align_orient];

      $totalRanges++;
    }
  }
  close IN;
  
  # Collapse redundancy and overlapping hits
  my @orientArr = ();
  my ( $tmpFH, $tmpFilename ) = tempfile( UNLINK => 1, SUFFIX => ".seqlist", DIR => "." );
  my $totalRegions = 0;
  foreach my $id ( keys %ranges ) {
    my $prevStart;
    my $prevEnd;
    my $prevOrient;
    print "ID = $id\n" if ( $DEBUG );
    # TODO: Track the count of overalapping ranges are on the +/- strand and use the max count to set the final orientation of the range
    foreach my $range ( sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$ranges{$id}} ) {
      my $start = $range->[0];
      my $end = $range->[1];
      my $orient = $range->[2];
      print "Reading: $start-$end $orient\n" if ( $DEBUG );
      if ( defined $prevStart )
      {
        next if ( $start == $prevStart && $end == $prevEnd );
        if ( $start <= $prevEnd + $gapTolerance ) {
          $start = $prevStart;
          print "overlapping elements\n" if ( $DEBUG );
        }else {
          push @orientArr, $prevOrient;
          # switch to zero-based, half-open
          print $tmpFH "$id:". ($prevStart-1) . "-" . $prevEnd . "\n";
          $totalRegions++;
          print "Saving: $prevStart-$prevEnd $prevOrient\n" if ( $DEBUG );
        }
      }
      $prevStart = $start;
      $prevEnd = $end;
      $prevOrient = $orient;
    }
    if ( defined $prevEnd )  {
      push @orientArr, $prevOrient;
      # switch to zero-based, half-open
      print $tmpFH "$id:". ($prevStart-1) . "-" . $prevEnd . "\n";
      $totalRegions++;
      print "Saving: $prevStart-$prevEnd $prevOrient\n" if ( $DEBUG );
    }
  }
  close $tmpFH;
  
  # Only report on what we found
  print "Multiple Alignment Edge Stats:\n";
  print "  Sequence Regions : $totalRegions\n";
  print "  Alignments* : $totalRanges\n";
  print "    Left  Edge: " . ( $totalRanges - $leftAdequate ) . " aligned sequence(s) needed extension, $leftUnextendable sequence(s) were unextendable\n";
  print "    Right Edge: " . ( $totalRanges - $rightAdequate ) . " aligned sequence(s) needed extension, $rightUnextendable sequence(s) were unextendable\n";
  print "  * There may be more than one alignment per sequence region\n";
  
  # Generate sequences
  my $cmd = "$TWOBITTOFA -seqList=$tmpFilename $genome stdout";
  open IN,"$cmd|" or die "Could not run $cmd!\n";
  if ( $outputFile ne "stdout" ) {
    open OUT,">$outputFile" or die "Could not open $outputFile for writing!\n";
  }
  my $idx = 0;
  my $seq = "";
  my $id;
  my $start;
  my $end;
  while (<IN>){
    if ( />(\S+)\:(\d+)-(\d+)/ ||
         />(\S+)/ ) {
      my $tmp_id = $1;
      my $tmp_start = $2;
      my $tmp_end = $3;
      $tmp_start = 0 if ( ! defined $tmp_start );
      $tmp_end = $seqLens{$tmp_id} if ( ! defined $tmp_end );
      if ( $seq ) {
        my $outID = $id . "_" . ($start + 1) . "_$end";
        if ( $orientArr[$idx] eq "-" ) {
          $seq = reverse( $seq );
          $seq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;
          $outID .= "_R";
        }
        if ( $outputFile eq "stdout" ) {
          print ">$outID\n$seq\n";
        }else {
          print OUT ">$outID\n$seq\n";
        }
        $idx++;
      }
      $id = $tmp_id;
      $start = $tmp_start;
      $end = $tmp_end;
      $seq = "";
      next;
    }
    s/[\n\r\s]+//g;
    $seq .= uc($_);
   }
   if ( $seq ) {
    my $outID = $id . "_" . ($start + 1) . "_$end";
    if ( $orientArr[$idx] eq "-" ) {
      $seq = reverse( $seq );
      $seq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;
      $outID .= "_R";
    }
    if ( $outputFile eq "stdout" ) {
      print ">$outID\n$seq\n";
    }else {
      print OUT ">$outID\n$seq\n";
    }
  }
  close IN;
  close OUT if ( $outputFile ne "stdout" );
  
  unless( $DEBUG ) {
    unlink( $tmpFilename ) if ( -e $tmpFilename );
  }
}


1;
