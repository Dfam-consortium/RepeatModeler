#!/usr/bin/perl -w

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
use File::Temp qw/ tempfile tempdir /;
#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use NCBIBlastSearchEngine;
use SearchResult;
use SearchResultCollection;


my $ucscToolsDir = "/usr/local/bin";
my $genome_file = "foo.2bit";


# Obtain seq lengths -- do this once
my %seqLens = ();
open IN,"$ucscToolsDir/twoBitInfo $genome_file stdout|" or die "Could not run $ucscToolsDir/twoBitInfo on $genome_file!\n";
my $lines = 0;
while ( <IN> ) {
  $lines++;
  if ( /^(\S+)\s+(\d+)/ ) {
    $seqLens{$1} = $2;
  }
}
close IN;


&addFlankingSequence($genome_file, "out", "newrepseq", \%seqLens, 100, 100);
# Do something
&addFlankingSequence($genome_file, "out", "newrepseq", \%seqLens, 100, 100);
# Do something



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
  my $cntTermLeft = 0;
  my $cntTermRight = 0;
  my $cntExtLeft = 0;
  my $cntExtRight = 0;
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
      my $orient = "+";
      $orient = "-" if ( $flds[8] eq "C" );
  
      my $id = "";
      my $start;
      my $end;
      # These have traditionally been in 1-based full-closed coordinates
      if ( $raw_id =~ /^(\S+)\_(\d+)\_(\d+)\_?(R?)$/ )
      {
        $id = $1;
        $start = $2;
        $end = $3;
      }else
      {
        die "I don't know how to parse this id: $raw_id\n";
      }
  
      if ( ! exists $seqLens{$id} ) {
        warn "Could not find $id in 2bit file!\n";
      }
      my $slen = $seqLens{$id};
  
      # Does the current sequence orientation still make sense or has
      # it reverted back to the forward strand?
      if ( $raw_id =~ /_R/ || $orient eq "-" && 
           !($raw_id =~ /_R/ && $orient eq "-" ) ) 
      {
        # Reversed left/right
        my $leftExtBP = $flankleft - $align_start - 1;
        if ( $leftExtBP > 0 ) {
          $end += $leftExtBP;
          if ( $end > $slen ) {
            $end = $slen;
            $cntTermLeft++;
          }else {
            $cntExtLeft++;
          }
        }
        my $rightExtBP = $flankright - $align_remain;
        if ( $rightExtBP > 0 ) {
          $start -= $rightExtBP;
          if ( $start < 1 ) {
            $start = 1;
            $cntTermRight++;
          }else {
            $cntExtRight++;
          }
        }
        push @{$ranges{$id}}, [$start, $end, "-"];
      }else {
        # Normal left/right
        my $leftExtBP = $flankleft - $align_start - 1;
        if ( $leftExtBP > 0 ) {
          $start -= $leftExtBP;
          if ( $start < 1 ) {
            $start = 1;
            $cntTermLeft++;
          }else {
            $cntExtLeft++;
          }
        }
        my $rightExtBP = $flankright - $align_remain;
        if ( $rightExtBP > 0 ) {
          $end += $rightExtBP;
          if ( $end > $slen ) {
            $end = $slen;
            $cntTermRight++;
          }else {
            $cntExtLeft++;
          }
        }
        push @{$ranges{$id}}, [$start, $end, "+"];
      }
      $totalRanges++;
    }
  }
  close IN;
  
  # Collapse redundancy and overlapping hits
  my @orient = ();
  my ( $tmpFH, $tmpFilename ) = tempfile( UNLINK => 0, SUFFIX => ".seqlist", DIR => "." );
  my $totalRegions = 0;
  foreach my $id ( keys %ranges ) {
    my $prevStart;
    my $prevEnd;
    my $prevOrient;
    print "ID = $id\n" if ( $DEBUG );
    foreach my $range ( sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$ranges{$id}} ) {
      my $start = $range->[0];
      my $end = $range->[1];
      my $orient = $range->[2];
      print "Reading: $start-$end $orient\n" if ( $DEBUG );
      if ( defined $prevStart )
      {
        next if ( $start == $prevStart && $end == $prevEnd );
        if ( $start <= $prevEnd ) {
          $start = $prevStart;
          print "overlapping elements\n" if ( $DEBUG );
        }else {
          push @orient, $prevOrient;
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
      push @orient, $prevOrient;
      # switch to zero-based, half-open
      print $tmpFH "$id:". ($prevStart-1) . "-" . $prevEnd . "\n";
      $totalRegions++;
      print "Saving: $prevStart-$prevEnd $prevOrient\n" if ( $DEBUG );
    }
  }
  close $tmpFH;
  
  if ( 0 ) 
  {
    # Only report on what we found
    print "Multiple Alignment Edge Stats:\n";
    print "  Sequence Regions : $totalRegions\n";
    print "  Alignments : $totalRanges\n";
    print "    Left  Edge: $cntExtLeft aligned sequences needed extension, $cntTermLeft sequences were unextendable\n";
    print "    Right Edge: $cntExtRight aligned sequences  needed extension, $cntTermRight sequences were unextendable\n";
    print "  * There may be 1 or more alignments per sequence region\n";
  }else {
    # Generate sequences
    my $cmd = "$ucscToolsDir/twoBitToFa -seqList=$tmpFilename $genome stdout";
    open IN,"$cmd|" or die "Could not run $cmd!\n";
    open OUT,">$outputFile" or die "Could not open $outputFile for writing!\n";
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
          if ( $orient[$idx] eq "-" ) {
            $seq = reverse( $seq );
            $seq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;
            $outID .= "_R";
          }
          print OUT ">$outID\n$seq\n";
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
      if ( $orient[$idx] eq "-" ) {
        $seq = reverse( $seq );
        $seq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;
        $outID .= "_R";
      }
      print OUT ">$outID\n$seq\n";
    }
    close IN;
    close OUT;
  }
  
  unless( $DEBUG ) {
    unlink( $tmpFilename ) if ( -e $tmpFilename );
  }
}


1;
