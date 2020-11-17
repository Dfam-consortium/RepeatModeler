#!/usr/bin/perl 
##---------------------------------------------------------------------------##
##  File:
##      @(#) extend-stk.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
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

extend-stk - Use RepeatAfterMe/ExtendAlign to automatically extend seed
             alignments.


=head1 SYNOPSIS

extend-stk.pl  [-debug|d]
               -assembly|a <*.2bit>
               -input|i <*.stk>
               -output|o <*.stk>

=head1 DESCRIPTION

Use the RepeatAfterMe ExtendAlign tool to try and extend seed alignment
that likely fragements of longer TE.  This is often the case with the
output of de novo TE discovery programs such as RepeatScout and RECON
the primary algorithms used in RepeatModeler.  

The options are:

=over 4

=item -assembly <*.2bit>

The genome assembly in UCSC 2bit format which contains the sequences
referenced in the input seed alignment.

=item -input <*.stk>

A stockholm file containing on or more seed alignments.

=item -output <*.stk>

The file to hold the results of the seed alignments after extension.

=back

=head1 ALSO

RepeatMasker, RepeatModeler, Dfam.org

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use RepModelConfig;
use MultAln;
use EMBL;
use SeedAlignment;
use SeedAlignmentCollection;
use File::Path 'rmtree';
use File::Temp qw/ tempfile tempdir /;
use Cwd;


# Program version
my $Version = 0.1;

#
# Paths
#
my $RepeatAfterMeDir = "/home/rh105648e/RepeatAfterMe";


# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-assembly|a=s',
                    '-input|i=s',
                    '-output|o=s',
                    '-debug|d'
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

unless ( $options{'assembly'} ) {
  print "\n\n-assembly is a required field!\n\n";
  usage();
}

unless ( $options{'output'} ) {
  print "\n\n-output is a required field!\n\n";
  usage();
}

unless ( $options{'input'} ) {
  print "\n\n-input is a required field!\n\n";
  usage();
}

my $genomeFile = $options{'assembly'};
my $inputFile = $options{'input'};
my $outputFile = $options{'output'};

# Open up the stockholm file and read in data
my $stockholmColl = SeedAlignmentCollection->new();
open my $IN, "<$inputFile"
     or die "Could not open up stockholm file $inputFile for reading!\n";
$stockholmColl->read_stockholm( $IN );
close $IN;

my $CLEANUP = 1;
$CLEANUP = 0 if ( $DEBUG );

open OUTPUT, ">$outputFile" or die;

my $failed = 0;
my $segmental = 0;
my $triggerSave = 0;
for ( my $i = 0 ; $i < $stockholmColl->size() ; $i++ ) {
  my $seedAlign   = $stockholmColl->get( $i );
  my $id          = $seedAlign->getId();
 
  # NOTE: This currently ignores the stockholm reference line and generates
  #       one by calling a consensus.
  my $tdiv = 14;
  my $mAlign = MultAln->new( seedAlignment => $seedAlign );
  my ( $null, $totDiv, $avgDiv ) = $mAlign->kimuraDivergence();
  print "MultAln Kimura divergence: $avgDiv\n";
  
  # TODO: make command line override
  $tdiv = $avgDiv;
  
  # This overrides the multaln method
  # Divergence is coming from generateSeedAlignment.pl:
  # Source:gsa, mDiv=23.22, GCA_000699945.1_ASM69994v1_genomic.fna.2bit:108
  if ( $seedAlign->getDescription() =~ /mDiv=(\d+\.\d+)/ )
  {
    $tdiv = $1;
  }

  my $div = 14;
  if ($tdiv >= 16) {
    $div = 18;
    if ($tdiv >= 19) {
      $div = 20;
      if ($tdiv >= 22.5) {
        $div = 25;
      }
    }
  }
 
  my $reference = $mAlign->getReferenceSeq();
  $reference =~ s/\-//g;

  my $tdir = tempdir( DIR => cwd, CLEANUP => $CLEANUP );
  print "Working on $id.. [ cons length = " . length($reference) . " ] in $tdir\n";
  
  open TSV,">$tdir/linup.tsv" or die "Could not open $tdir/linup.tsv for writing!\n";

  for ( my $j = 0; $j < $seedAlign->alignmentCount(); $j++ ){
    my ( $assemblyName, $sequenceName, $start, $end, $orient, $sequence ) =
      $seedAlign->getAlignment( $j );

    $start--; # convert to zero-based half open

    # insertions...
    if ( 0 ) {
      if ( $assemblyName ne "" ){
        $sequenceName = $assemblyName.":".$sequenceName;
      }
    }
    # If insertions against assembly...insertions against insertions is already fine
    if ( 0 ) {
      print "Warning accounting for 50bp padding for 200Mammals insertions\n";
      $start -= 50;
      $end -= 50;
    }

    # Tolerate up to 10bp from edge
    if ( $sequence =~ /^[\.]{0,10}[^\.]/ && $sequence =~ /[^\.][\.]{0,10}$/ ) {
      print TSV "$sequenceName\t".$start."\t$end\t1\t1\t$orient\n";
    }elsif ( $sequence =~ /^[\.]{0,10}[^\.]/ ) 
    {
      # can extend left
      print TSV "$sequenceName\t".$start."\t$end\t1\t0\t$orient\n";
    }elsif ( $sequence =~ /[^\.][\.]{0,10}$/ ) 
    {
      # can extend right
      print TSV "$sequenceName\t".$start."\t$end\t0\t1\t$orient\n";
    }else {
      print TSV "$sequenceName\t".$start."\t$end\t0\t0\t$orient\n";
    }
  } 
  close TSV;
  my $cmd = "$RepeatAfterMeDir/ExtendAlign -twobit $genomeFile -bandwidth 500 -matrix $div" . "p43g -ranges $tdir/linup.tsv -outtsv $tdir/repam-ranges.tsv -outfa $tdir/repam-repseq.fa -cons $tdir/repam-cons.fa >& $tdir/repam.log";
  print "Running extendalign..\n";
  print "Running: $cmd\n";
  system($cmd);
    
  #Extended right: 0 bp
  #Extended left : 0 bp
  #Program duration is 4.0 sec = 0.1 min = 0.0 hr
  my $t1 = `fgrep "Extended right: 0 bp" $tdir/repam.log`;
  my $t2 = `fgrep "Extended left : 0 bp" $tdir/repam.log`;
  my $newCons = "";
  #print "t1=>$t1< t2=>$t2<\n";
  if ( $t1 && $t2 ) {
    print "**Could not extend**\n";
    print OUTPUT "" . $seedAlign->toString();
    next;
  }elsif ( -e "$tdir/repam-cons.fa") {
    open IN,"<$tdir/repam-cons.fa" or die;
    my $ttid;
    my %seqs = ();
    while (<IN>){
      if ( />(\S+)/ )
      {
        $ttid = $1;
        next;
      }
      s/[\n\r\s]+//g;
      $seqs{$ttid} .= $_;
    }
    close IN;
    $newCons = $seqs{'left-extension'} . $reference . $seqs{'right-extension'};
    #print "newCons = $newCons\n";
  }else {
    print "Something went wrong with extend align:\n";
    print "  $cmd\n";
    $failed++;
    $triggerSave = 1;
    next;
  } 

  ##
  ## 
  ##
  $t1 = `fgrep "Extended right:" $tdir/repam.log`;
  $t2 = `fgrep "Extended left :" $tdir/repam.log`;
  my $rightExt = 0;
  if ( $t1 =~ /Extended\s+(left|right)\s*:\s+(\d+)/ ){
    $rightExt = $2;
  }
  my $leftExt = 0;
  if ( $t2 =~ /Extended\s+(left|right)\s*:\s+(\d+)/ ){
    $leftExt = $2;
  }
  print "Captured Extensions: left $leftExt bp, right $rightExt bp\n";
  if ( $leftExt > 9990 && $rightExt > 9990 )
  {
    print "Extension hit limits in both directions...probably a segmental duplication...keeping unextended\n";
    $segmental++;
    my $desc = $seedAlign->getDescription();
    $desc =~ s/[\n\r]+//g;
    $desc .= " [possibly part of segmental duplication]";
    $seedAlign->setDescription( $desc );
    print OUTPUT "" . $seedAlign->toString();
    next;
  }

  open OUT,">$tdir/repam-newrep.fa" or die;
  print OUT ">repam-newrep\n$newCons\n";
  close OUT;
  print "  extended by " . (length($newCons)-length($reference)) . " bp\n";
  my $prevCons;

  print "Refining consensus...\n";
  my $maxIter = 10;
  system("$FindBin::Bin/alignAndCallConsensus.pl -c $tdir/repam-newrep.fa -e $tdir/repam-repseq.fa -refine $maxIter -st -q");

  # Now rebuild the stockholm entry and save to ...
  my $singleSeedColl = SeedAlignmentCollection->new();
  open my $INN, "<$tdir/repam-newrep.stk"
    or die "Could not open up stockholm file $tdir/repam-newrep.stk for reading!\n";
  $singleSeedColl->read_stockholm( $INN );
  close $INN;
  my $finalAlign   = $singleSeedColl->get( 0 );
  $finalAlign->setId( $seedAlign->getId() );
  $finalAlign->setDescription( $seedAlign->getDescription() );
  $finalAlign->clearCitations();
  # Correct a strange mapping artefact in RepeatModeler data
  my $tClass = $seedAlign->getClassification();
  if ( $tClass eq "Unknown" ) {
    $tClass = "Interspersed_Repeat;Unknown";
  }
  $finalAlign->setClassification( $tClass );
  $finalAlign->clearClades();
  for ( my $k = 0; $k < $seedAlign->cladeCount(); $k++ ){
    $finalAlign->addClade( $seedAlign->getClade($k) );
  }
  #open OSTK,">$id.stk-final" or die;
  #print OSTK "" . $finalAlign->toString();
  #close OSTK;

  print OUTPUT "" . $finalAlign->toString();

  # Cleanup
  unless($DEBUG) {
    rmtree([ $tdir ]);
  }
}
close OUTPUT;


sub readConsensusFromFile {
   my $file = shift;
   open IN,"<$file" or die;
   my $seq;
   while (<IN>){
     if ( />\S+/ )
     {
       next;
     }
     s/[\n\r\s]+//g;
     $seq .= $_;
   }
   close IN;
   return $seq;
}


1;
