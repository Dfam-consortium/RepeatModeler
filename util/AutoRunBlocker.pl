#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) AutoRunBlocker.pl
##  Author:
##      Arian Smit         asmit@systemsbiology.org
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Automatically identify insertions/deletions in a Linup consensus
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
#******************************************************************************

=head1 NAME

AutoRunBlocker.pl - Automatically identify insertions/deletions in a Linup consensus

=head1 SYNOPSIS

  AutoRunBlocker.pl -linup|-l <ali file>
                -windowSize|-w #
                -minCopyAgreement|-mc #
                [-prompt|-p]
                [-minRatioAgreement|-mr #]

  or 

  AutoRunBlocker.pl -l <linup_ali_file> -w <windowSize> -mc <minCopyAgreement> 
 
  or

  AutoRunBlocker.pl -prompt -l <linup_ali_file> -w <windowSize> -mc <minCopyAgreement> -mr <minRatioAgreement>

=head1 DESCRIPTION

Runs Blocker.pl on an Linup format alignment file for every position
with a given window.  The results are then filtered using the provided
copy/ratio/N criteria and overlapping results are merged into distinct
unstable regions.  These larger regions are then fed to Blocker again
to identify the new consensus for the entire region.  If these larger
regions also meet the filter criteria the consensus is updated.

The options are:

=over 4

=item -linup <ali file>

A alignment file produced by Linup ( in "ali" standard format ).

=item -prompt 

If this option is used, each clustered result is displayed and the user is
prompted to either include it or exclude it in the consensus update.

=item -minRatioAgreement #

The minimum threshold of the ratio of count of majority length subsequences to
the number retaining the original length.  Default: 1

=item -windowSize #

The size of the window in bp.

=item -minCopyAgreement #

The minimum number of sequences that must agree on length over the window.

=back

=head1 SEE ALSO

ReapeatModeler, alignAndCallConsensus.pl etc.. 

=head1 COPYRIGHT

Copyright 2019-2021 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Arian Smit <asmit@systemsbiology.org>

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use Getopt::Long;
use POSIX qw(:sys_wait_h ceil floor);
use File::Copy;
use File::Spec;
use File::Path;
use File::Basename;
use Cwd qw(abs_path getcwd cwd);
use Data::Dumper;
#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};

# Program version
my $Version = $RepModelConfig::VERSION;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-linup|l=s',
                    '-windowSize|w=i',
                    '-prompt|p',
                    '-debug|d',
                    '-minCopyAgreement|mc=i',
                    '-minRatioAgreement|mr=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit( 1 );
}

my $DEBUG = 0;
$DEBUG = 1 if ( $options{'debug'} );

my $alifile;
my $window;
my $copymin;
my $ratio = 1;
if ( $options{'linup'} && $options{'windowSize'} && $options{'minCopyAgreement'} ) {
  $alifile = $options{'linup'};
  $window = $options{'windowSize'};
  $copymin = $options{'minCopyAgreement'};
  $ratio = $options{'minRatioAgreement'} if ( $options{'minRatioAgreement'} );
}elsif ( @ARGV >= 3 ) {
  $alifile = shift @ARGV;
  $window = shift @ARGV;
  $copymin = shift @ARGV;
  $ratio = shift @ARGV if ( @ARGV );
}else {
  usage();
}

if ( ! defined $window ) {
  print "\n\nWindow size is a required parameter!\n\n";
  usage;
}
if ( ! defined $copymin ) {
  print "\n\nMin copy agreement is a required parameter!\n\n";
  usage;
}


open (IN, $alifile) or die "Could not open up $alifile for reading!";
my @starts;
while (<IN>) {
  if (/consensus\s+(\d+)/ ) {
    push @starts, $1;
  }
}
close IN;
my %blocks = ();
foreach my $start (@starts) {
  my $name = "block$start";
  $blocks{$name} = [];
}
open (IN, $alifile) or die;
my $consensus = "";
while (<IN>) {
  if (/consensus\s+(\d+)\s+(\S+)/ ) {
    $consensus .= $2;
    my $name = "block$1";
    my @string = split "", $2;
    for (my $i = 0; $i <= $#string; ++$i) {
      if ($string[$i] =~ /\w/) {
        # Columns which contain a consensus base and not a gap
        push @{$blocks{$name}}, $i; #still zero-based
      }
    }
  }
}
$consensus =~ s/-//g;

# Gather blocker output
my @blockerData = ();
while ($starts[0]) {
    my $start = shift @starts;
    my $name = "block$start";
    # going through all positions in the alignment where the consensus has a base
    for (my $j = 0; $j <= $#{$blocks{$name}}; ++$j) {
      my $k = $j + $window - 1;
      my $beg = $blocks{$name}->[$j] + 1; # blocker.pl is one based
      my $begin = $start + $j;
      my $eind = $begin + $window - 1;

      if ($blocks{$name}->[$k]) {  # the end position is still found in the same block
        my $end = $blocks{$name}->[$k] + 1;
        my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
           &callBlocker( $alifile, $start, $beg, $start, $end);
        #print "Calling Same Block [ $start, $beg, $start, $end ]\n";
        if ( $numMajorityLen ) {
          print "Blocker: $begin-$eind [$start,$beg,$start,$end]: $numMajorityLen, $numOrigLen, $newseqlen, $ncount : $newseq\n" if ( $DEBUG );
          push @blockerData, [$numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount, $begin, $eind, $start, $beg, $start, $end];
        }
      } else { # the end positions is in a next block
        my @startsleft = @starts; # current start already shifted from @starts
        my $remain = $k - $#{$blocks{$name}} - 1; 
        #print "Remain = $remain  k= $k j=$j window=$window number of positions in block = " . $#{$blocks{$name}} . "\n";
        while ($remain >= 0 && $startsleft[0]) {
          my $nxtstart = shift @startsleft;
          my $naam = "block$nxtstart";
          if (defined $blocks{$naam}->[$remain]) { #can be 0
            my $end = $blocks{$naam}->[$remain] + 1; 
            #print "Calling [ $start, $beg, $nxtstart, $end ]\n";
            my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount) =
               &callBlocker( $alifile, $start, $beg, $nxtstart, $end );
            if ( $numMajorityLen ) {
              print "Blocker: $begin-$eind [$start,$beg,$nxtstart,$end]: $numMajorityLen, $numOrigLen, $newseqlen, $ncount : $newseq\n" if ( $DEBUG );
              push @blockerData, [$numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount, $begin, $eind, $start, $beg, $nxtstart, $end];
            }
            $remain = -1;
          } else {
            $remain -= ($#{$blocks{$naam}} + 1);
          }
        } # end while
      }
    }
  }


my $prevRec;
my $allowedGapDist = $window/5; #
my $clusterStartBlock;
my $clusterStartColumn;
my $clusterStartCons;
my @changes = ();
foreach my $bRec ( @blockerData ) {
  # Records: [$numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount, $consbegin, $consend, $startBlock, $startColumn, $endBlock, $endColumn]
  if ( ! defined $clusterStartBlock ) {
    $clusterStartBlock = $bRec->[7];
    $clusterStartColumn = $bRec->[8];
    $clusterStartCons = $bRec->[5];
  }

  if ( $prevRec ) {
    if ( $bRec->[5] - $allowedGapDist - 1 <= $prevRec->[6] ) {
      # Clustered with previous
    }else {
      # Independent of previous
      if ( defined $clusterStartBlock ) {
        # Emit cluster
        print "Cluster: ./Blocker ali $clusterStartBlock $clusterStartColumn " . $prevRec->[9] . " " . $prevRec->[10] . "\n" if ( $DEBUG );
        my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount, $blockerout) =
             &callBlocker( $alifile, $clusterStartBlock, $clusterStartColumn, $prevRec->[9], $prevRec->[10]);
        if ( $numMajorityLen ) {
          if ( $options{'prompt'} ) {
            print STDERR "-----------------------------------------------------------\n";
            print STDERR "Blocker Results from pos $clusterStartCons:\n";
            print STDERR "$blockerout\n";
            print STDERR "Accept (Y/N)? [Y]: "; 
            my $answer = <STDIN>;
            $answer =~ s/[\n\r]+//g;
            if ( $answer =~ /y|yes/i || $answer eq "" ) {
              push @changes, [$clusterStartCons, $prevRec->[6], $newseq];
            }
          }else {
            print "Success: $clusterStartCons-$prevRec->[6]:$newseq\n" if ( $DEBUG );
            push @changes, [$clusterStartCons, $prevRec->[6], $newseq];
          }
        }
      }
      $clusterStartBlock = $bRec->[7];
      $clusterStartColumn = $bRec->[8];
      $clusterStartCons = $bRec->[5];
    }
  }
  $prevRec = $bRec; 
}
if ( defined $clusterStartBlock ) {
  print "Cluster: ./Blocker ali $clusterStartBlock $clusterStartColumn " . $prevRec->[9] . " " . $prevRec->[10] . "\n" if ( $DEBUG );
  my ( $numMajorityLen, $numOrigLen, $newseq, $newseqlen, $ncount, $blockerout) =
     &callBlocker( $alifile, $clusterStartBlock, $clusterStartColumn, $prevRec->[9], $prevRec->[10]);
  if ( $numMajorityLen ) {
      if ( $options{'prompt'} ) {
        print STDERR "-----------------------------------------------------------\n";
        print STDERR "Blocker Results from pos $clusterStartCons:\n";
        print STDERR "$blockerout\n";
        print STDERR "Accept (Y/N)? [Y]: "; 
        my $answer = <STDIN>;
        $answer =~ s/[\n\r]+//g;
        if ( $answer =~ /y|yes/i || $answer eq "" ) {
          push @changes, [$clusterStartCons, $prevRec->[6], $newseq];
        }
      }else {
        print "Success: $clusterStartCons-$prevRec->[6]:$newseq\n" if ( $DEBUG );
        push @changes, [$clusterStartCons, $prevRec->[6], $newseq];
      }
  }
}

if ( @changes ) 
{
  print STDERR "There are " . scalar(@changes) . " changes\n";
}else {
  print STDERR "There are no changes\n";
}

my $newconsensus = $consensus; 
foreach my $change ( sort { $b->[0] <=> $a->[0] } @changes ) {
  substr($newconsensus, $change->[0] - 1, $change->[1] - $change->[0] + 1) = $change->[2]; 
}
print "Consensus:\n$consensus\n" if ( $DEBUG );
print "New Consensus:\n$newconsensus\n" if ( $DEBUG );

print ">rep\n$newconsensus\n";




sub callBlocker {
  my $linup = shift;
  my $startBlock = shift;
  my $startColumn = shift;
  my $endBlock = shift;
  my $endColumn = shift;

  my $blockerout;
  if ( ! defined $endColumn ) {
   $blockerout = `$FindBin::RealBin/Blocker.pl $linup $startBlock $startColumn $endBlock`;
  }else {
   $blockerout = `$FindBin::RealBin/Blocker.pl $linup $startBlock $startColumn $endBlock $endColumn`;
  }

  my $numMajorityLen = 0;
  my $numOrigLen = 0;
  my $newseq = "";
  my $newseqlen = 0;
  my $ncount = 0;
  # Length difference. 66 of 113 for 11 bp (33 second best for 10 bp, 33 for original 10 bp)
  if ( $blockerout =~ /difference\.\s*(\d+).+\s(\d+)\s+for\s+original/ ) {
     $numMajorityLen = $1;
     $numOrigLen = $2;
  }
  if ( $blockerout =~ /new\s+([ACGTN]+)/ ) {
     $newseq = $1;
     $newseq =~ s/^N+//; #otherwise the beg/end are not reported
     $newseq =~ s/N+$//;
     $ncount = ($newseq =~ tr/N//);
  }
  if ( $blockerout =~ /Sequence difference[\.\s]+(\d+)/) {
     $numMajorityLen = $1;
  }
  my $newlength = length($newseq);
  my $adjustedratio = $ratio;
  my $nratiocut = 10;
  if ($newlength > $window) {
    $adjustedratio -= 0.5;
    $nratiocut -= ($newlength - $window)/4;
  }
  if ( $numMajorityLen >= $copymin && (!$numOrigLen || $numMajorityLen/$numOrigLen >= $adjustedratio) 
       && ($ncount < 2 || $newlength/$ncount > $nratiocut ) ) 
  {
    return ( $numMajorityLen, $numOrigLen, $newseq, $newlength, $ncount, $blockerout);
  }else {
    return ( 0, 0, "", 0, 0);
  }
}


