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

 extendFlankingSeqs.pl [-version] [[-r(ight_flank) #] [-l(eft_flank) #] | [-f(lank) #]]
                       [-g(ap) #] [-r(eindex)]
                       -d(atabase) <2bit file> 
                       -i(nput) <cross_match file> -o(utput) <fasta file>

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

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
                    '-reindex',
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
  print "\n\nDatabase doesn't exist!  Must supply a TwoBit file!\n\n";
  usage();
}

if ( ! exists $options{'input'} || ! -s $options{'input'} ){
  print "\n\nMust supply input alignments in crossmatch format!\n\n";
  usage();
}

if ( ! exists $options{'output'} ){
  print "\n\nMust supply output filename!\n\n";
  usage();
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

if ( ! -e $indexFile || $options{'reindex'} ) {
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
      if ( $raw_id =~ /^(\S+)\_(\d+)\_(\d+)\_?(R?)$/ )
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
        }elsif ( $4 eq "R" ) {
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

print "SeqRecord: $seqRecordID:$seqRecordStart-$seqRecordEnd ($seqRecordOrient) <=> Alignment: $align_start-$align_end [$align_remain] ($align_orient)\n";
  
      $leftAdequate++ if ( $align_start-1 >= $flankleft );
      $rightAdequate++ if ( $align_remain >= $flankright ); 

      my $globalStart;
      my $globalEnd;
      if ( $seqRecordOrient ne $align_orient ) {
        if ( $seqRecordOrient eq "-" ) {
          #  rec:200-100  hit:1-10 => 190-200 (-)
          $globalStart = $seqRecordEnd - $align_end + 1; # right
          $globalEnd = $seqRecordEnd - $align_start + 1; # left
        }else {
          #  rec:100-200 hit:10-1 => 100-110 (-)
          $globalStart = $seqRecordStart + $align_start - 1; # right
          $globalEnd = $seqRecordStart + $align_end - 1; # left
        }
        $rightUnextendable++ if ( $globalStart == 1 );  
        if ( $align_remain < $flankright ) {
          $globalStart -= $flankright;
          $globalStart = 1 if ( $globalStart < 1 );
        }
        $leftUnextendable++ if ( $globalEnd == $seqRecordLen );  
        if ( $align_start-1 < $flankleft ) {
          $globalEnd += $flankleft;
          $globalEnd = $seqRecordLen if ( $globalEnd > $seqRecordLen );
        }
print " Range: $globalStart-$globalEnd (-)\n";
        push @{$ranges{$seqRecordID}}, [$globalStart, $globalEnd, "-"];
      }else {
         if ( $seqRecordOrient eq "-" ) {
          #  rec:200-100  hit:10-1 => 190-200 (+)
          $globalStart = $seqRecordEnd - $align_end + 1; # left
          $globalEnd = $seqRecordEnd - $align_start + 1; # right
print "inter: $globalStart-$globalEnd $seqRecordEnd\n";
        }else {
          #  rec:100-200 hit:1-10 => 100-110 (+)
          $globalStart = $seqRecordStart + $align_start - 1; # left
          $globalEnd = $seqRecordStart + $align_end - 1; # right
        }
        $leftUnextendable++ if ( $globalStart == 1 );  
        if ( $align_start-1 < $flankleft ) {
          $globalStart -= $flankleft;
          $globalStart = 1 if ( $globalStart < 1 );
        }
        $rightUnextendable++ if ( $globalEnd == $seqRecordLen );  
        if ( $align_remain < $flankright ) {
          $globalEnd += $flankright;
          $globalEnd = $seqRecordLen if ( $globalEnd > $seqRecordLen );
        }
print " Range: $globalStart-$globalEnd (+)\n";
        push @{$ranges{$seqRecordID}}, [$globalStart, $globalEnd, "+"];
      }
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
  print "  Alignments : $totalRanges\n";
  print "    Left  Edge: " . ( $totalRanges - $leftAdequate ) . " aligned sequences needed extension, $leftUnextendable sequences were unextendable\n";
  print "    Right Edge: " . ( $totalRanges - $rightAdequate ) . " aligned sequences needed extension, $rightUnextendable sequences were unextendable\n";
  print "  * There may be 1 or more alignments per sequence region\n";
  if ( 1 ) { 
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
  }
  
  unless( $DEBUG ) {
    unlink( $tmpFilename ) if ( -e $tmpFilename );
  }
}


1;
