#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) alignAndCallCons2.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Arian Smit         asmit@systemsbiology.org
##  Description:
##      A script to generate a MSA from a set of TE copies and call
##      a new consensus.
##

=head1 NAME

alignAndCallCons2.pl - Align TE instances and call consensus

=head1 SYNOPSIS

  alignAndCallCons2.pl -e(lements) <*.fa> | -de(faults) 
                      [-help]

  Example:

    ./alignAndCallCons2.pl -de

        Assumes that a file named "repseq" exists and contains
        all the sequences to be aligned.  This is a deprecated option
        will eventually be removed in favor of the -e parameter.

  or

    ./alignAndCallCons2.pl -e myinstances.fa

=head1 DESCRIPTION

The options are:

=over 4

=item -e(lements) <*.fa>

A multi-sequence FASTA file containing sequences believed to be elements of 
a single TE family.  

=back

=head1 SEE ALSO

ReapeatModeler, dothemsimple.pl, extendcons.pl, createxmoutstk.pl etc..

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
use Data::Dumper;
use Getopt::Long;
use Cwd;
use File::Spec;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use NCBIBlastSearchEngine;
use CrossmatchSearchEngine;
use SearchResult;
use SearchResultCollection;

# Program version
my $Version = 0.1;

#
# Paths
#
my $phrapDir = "/usr/local/phrap";
my $matrixDir = "$FindBin::RealBin/../Matrices";
my $RMBLAST_DIR = $RepModelConfig::configuration->{'RMBLAST_DIR'}->{'value'};
my $defaultEngine = "crossmatch";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-consensus|c=s',
                    '-elements|e=s',
                    '-divergencemax=s', 
                    '-defaults|de',
                    '-finishedext|f=s',
                    '-matrix|ma=s',
                    '-minmatch|mi=s',
                    '-score|sc=s',
                    '-stockholm|st',
                    '-outdir=s',
                    '-include_reference|inc',
                    '-interactive|int',
                    '-prune|p=s',
                    '-quoteprune|qu=s',
                    '-html|h',
                    '-quiet|q',
                    '-refine|re:s',
                    '-bandwidth|b=s',
                    '-crossmatch|cr',
                    '-rmblast|rm',
                    '-help',
                    '-filter_overlapping|fi',
                    '-onlyBestAlignment', # Deprecated, now 'filter_overlapping'
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

usage() unless ( ! $options{'help'} && (
                 $options{'defaults'} || ( $options{'consensus'} && $options{'elements'} )));


# Output directory
my $outdir = cwd;

#
# Process the consensus file:
#    - Unless specifified default to "rep"
#
#    - If it contains more than one sequence we assume that we are 
#      in a mode to compete and build several subfamilies against each
#      other ( aka. dothemmultiple ).
#
#    - If a consensus sequence contains the suffix "#buffer" it is treated
#      as a decoy for attracting instances and will not be refined.
#      
my $conFile = "rep";
if ( exists $options{'consensus'} ) {
  if ( -s $options{'consensus'} ) {
    $conFile = $options{'consensus'}; 
    my($filename, $dirs, $suffix) = fileparse($options{'consensus'});
    $outdir = $dirs;
  }else {
    print "\n\nConsensus file \"$options{'consensus'}\" is missing or empty! Use '-help' for option help.\n\n";
    exit(1);
  }
}elsif ( ! -s $conFile )
{
  print "\n\nDefault consensus file \"rep\" is missing or empty! See -c option for details. Use '-help' for option help.\n\n";
  exit(1);
}
my ( $numRefineableCons, $numBuffers, $consRecs ) = &processConFile( $conFile );
if ( $numRefineableCons == 0 && $numBuffers > 0 ) {
  print "\n\nThe consensus file only contains buffer sequences.  There must be at least one non-buffer sequence. Use '-help' for option help.\n\n";
  exit(1);
}

# Check to see if there are already some previous runs in this directory.
my $firstFreeBackupIdx = 0;
for ( my $i = 1; $i < 25; $i++ ) {
  if ( ! -e "$conFile.$i" ) {
    $firstFreeBackupIdx = $i;
    last;
  }
}


# Override all other methods at picking an outdir if this option is provided.
$outdir = $options{'outdir'} if ( $options{'outdir'} );
# Assumption below is that the directory doesn't already include trailing "/"
$outdir =~ s/\/$//g;
if ( ! -d $outdir ) {
  die "\n\nOutput directory $outdir doesn't exist!\n\n";
}

my $elesFile = "repseq";
if ( exists $options{'elements'} ) {
  if ( -s $options{'elements'} ) {
    $elesFile = $options{'elements'}; 
  }else {
    print "\n\nElements file \"$options{'elements'}\" is missing or empty!\n\n";
    exit(1);
  }
}elsif ( ! -s $elesFile )
{
  print "\n\nDefault elements file \"repseq\" is missing or empty! See -e option for details.\n\n";
  exit(1);
}

my $engine = $defaultEngine;
$engine = "rmblast" if ( $options{'rmblast'} );
$engine = "crossmatch" if ( $options{'crossmatch'} );
my $matrixPath;
if ( $engine eq "rmblast" ) {
  $matrixPath = "$matrixDir/ncbi/nt";
}else{
  $matrixPath = "$matrixDir/crossmatch";
}

# Configure the matrix
my $matDiv = "25";
my $matGC = "41";
if ($options{'matrix'}) {
  # Simple single integer form: e.g. -matrix 25
  if ($options{'matrix'} =~ /^(\d+)$/) {
    $matDiv = $1;
  }elsif ($options{'matrix'} =~ /^(\d{2})p(\d{2})g?$/) {
    $matDiv = $1;
    $matGC = $2;
  }else{
    die "Cannot parse matrix parameter: $options{'matrix'}  Expected either a single integer or '##p##g' format.\n";
  }
}

my $matSpec;
my $gapInit;
my $insGapExt;
my $delGapExt;
if ( $matDiv =~ /14|18|20|25/ ){
  if ( $matDiv == 25 ){
    $gapInit = -25; $insGapExt = -5; $delGapExt = -4;
  }elsif ( $matDiv == 20 ) {
    $gapInit = -28; $insGapExt = -6; $delGapExt = -5;
  }elsif ( $matDiv == 14 ) {
    $gapInit = -33; $insGapExt = -7; $delGapExt = -6;
  }else {
    $gapInit = -30; $insGapExt = -6; $delGapExt = -5;
  }
  if ( $matGC =~ /37|39|41|43|45|49|51|53/ ) {
    $matSpec = $matDiv . "p" . $matGC . "g-Hpad.matrix";
  }else {
    die "Error: Matrix parameter can only choose from 37, 39, 41, 43, 45, 49, 51 or 53 for matrix GC levels\n";
  }
}else{
  die "Error: Matrix parameter can only choose from 14, 18, 20 or 25 for matrix divergence levels\n";
}
if ( ! -s "$matrixPath/$matSpec" ) {
  die "Error: Matrix parameter resolved to $matrixPath/$matSpec, which doesn't exist!\n";
}

my $minmatch  = 7;
$minmatch = $options{'minmatch'} if ( exists $options{'minmatch'} );

my $minscore = 200;
$minscore = $options{'score'} if ( exists $options{'score'} );

my $bandwidth = 40;
$bandwidth = $options{'bandwidth'} if ( exists $options{'bandwidth'} );

my $maxdiv = 60;
$maxdiv = $options{'divergencemax'} if ( exists $options{'divergencemax'} );

my $maxRefineIterations = 5;
$maxRefineIterations = $1 if ( exists $options{'refine'} && $options{'refine'} =~ /(\d+)/ );

my ($ext5done, $ext3done) = ();
if ($options{'finishedext'}) {
  if ($options{'finishedext'} =~ /^5/) {
    $ext5done = 1;
  } elsif ($options{'finishedext'} =~ /^3/) {
    $ext3done = 1;
  } elsif($options{'finishedext'} =~ /^[bB]/) {
    $ext3done = 1;
    $ext5done = 1;
  } else {
    die "-f $options{'finishedext'} not recognized\n";
  }
}

my $prunecutoff = 0;
my $pruneratio = 0;
my $npad = 0;
# the alignment of sequences with Npads at query termini and Hpads at
# consensus termini often involves Ns aligning to the real consensi termini
# due to the attraction of the Hpads (N vs H scores 7 like anything else)
# these should not count to the number of people aligned at those positions
if ($options{'prune'}) {
  $prunecutoff = $options{'prune'};
  if ($prunecutoff =~ s/[nN](\d+)$//) {
    $npad = $1;
  }
}
if ($options{'quoteprune'}) {
  die "quoteprune is NOT implemented...\n";
  $pruneratio = $options{'quoteprune'};
  if ($pruneratio =~ s/[nN](\d+)$//) {
    $npad = $1 unless $npad;
  }
}


##
## Main
##
unless ( $options{'quiet'} ) {
  print "##\n";
  print "## alignAndCallConsensus (aka dothemsimple/dothemmultiple.pl)\n";
  print "## Version $Version\n";
  print "##\n";
  if ( $numRefineableCons > 1 ) {
    print "# Multiple Subfamily Mode ";
  }else {
    print "# Single Family Mode ";
  }
  if ( $numBuffers > 0 ) {
    print "[ $numBuffers buffer sequence(s) ]";
  }
  print "\n";
  print "# Engine: $engine  Matrix: $matSpec, Bandwidth: $bandwidth,\n";
  print "#                  Minmatch: $minmatch, Minscore: $minscore,\n";
  print "#                  Maxdiv: $maxdiv, GapInit: $gapInit,\n";
  print "#                  InsGapExt: $insGapExt, DelGapExt: $delGapExt\n";
  if ( $options{'refine'} ) {
    print "# Refine: max $maxRefineIterations iterations\n";
  }else {
    print "# Refine: false\n";
  }
  if ( $firstFreeBackupIdx > 1 ) {
    print "# Starting Round Index: $firstFreeBackupIdx\n";
  }
}

my $searchEngineN;
if ( $engine eq "crossmatch" ) {
  $searchEngineN = CrossmatchSearchEngine->new(
                         pathToEngine=>"$phrapDir/cross_match" );
  $searchEngineN->setBandwidth( $bandwidth );
}else {
 $searchEngineN = NCBIBlastSearchEngine->new(
                       pathToEngine => "$RMBLAST_DIR/rmblastn" );
  # This parameter is a work-in-progress for rmblastn.  See the NCBIBlastSearchEngine for
  # details.
  # $searchEngineN->setBandwidth( -50 );
  $searchEngineN->setCores( 10 );
}

$searchEngineN->setGenerateAlignments( 1 );
$searchEngineN->setMatrix("$matrixPath/$matSpec");
# Default for crossmatch
$searchEngineN->setMaskLevel( 80 );
$searchEngineN->setMinScore( $minscore );
$searchEngineN->setGapInit( $gapInit );
$searchEngineN->setInsGapExt( $insGapExt );
$searchEngineN->setDelGapExt( $delGapExt );
$searchEngineN->setMinMatch( $minmatch );
$searchEngineN->setScoreMode( SearchEngineI::complexityAdjustedScoreMode );
$searchEngineN->setQuery( $elesFile );
# Consider making a NCBIBlastSearchEngine.pm method for using the "-sequence" parameter
# when only one sequence is in the database -- as is the case here.
$searchEngineN->setSubject( $conFile );

my $iterations = 0;
my $backupIdx = 1;
while ( 1 ) { 

  unless ( $engine eq "crossmatch" ) {
    # Always make a database as the sequence may have changed
    system(   "$RMBLAST_DIR/makeblastdb -blastdb_version 4 -out $conFile "
            . "-parse_seqids -dbtype nucl -in $conFile >& /dev/null ");
  }
 
  my %collected = ();
  #print STDERR "Search Parameters: " . $searchEngineN->getParameters() . "\n";
  my ( $status, $resultCollection ) = $searchEngineN->search();
  if ( $status )
  {
    print STDERR "Search Parameters: " . $searchEngineN->getParameters() . "\n";
    die "\nERROR from search engine (", $? >> 8, ") \n";
  } 
  if ( $engine eq "crossmatch" && -e "$elesFile.log" ) {
    unlink("$elesFile.log");
  }

  # Sort the results collection
  $resultCollection->sort(
    sub ($$) {
      $_[ 1 ]->getScore() <=> $_[ 0 ]->getScore();
     }
  );

  #
  # pruning
  #
  if ( $prunecutoff || $pruneratio ) {
  
    # Go through collection and determine how much to clip
    # consensi.  NOTE: This counts multiple alignments to
    # a single sequence as one contribution from the min to
    # max consensus range -- this takes care of tandem repeats
    # that would otherwise contribute too many counts.
    my %occBegin = ();
    my %occEnd = ();
    for ( my $k = 0 ; $k < $resultCollection->size() ; $k++ ) {
      my $resultRef = $resultCollection->get( $k );
      my $sName = $consRecs->{$resultRef->getSubjName()};
      $sName = $1 if ( $sName =~ /(\S+)\#.*/ );
      next if ( $consRecs->{$sName}->{'buffer'} == 1 );
      my $qName = $resultRef->getQueryName();
      my $qStart = $resultRef->getQueryStart();
      my $qRem = $resultRef->getQueryRemaining();
      my $min5 = 0;
      my $min3 = 0;
      if ( $npad ) {
        if ( $qStart <= $npad ) {
          $min5 = 1 + $npad + $qStart;
        }
        if ( $qRem < $npad ) {
          $min3 = $npad - $qRem;
        }
      }
      my $tbegin = $resultRef->getSubjStart() + $min5;
      my $tend = $resultRef->getSubjEnd() - $min3;
      if ( $resultRef->getOrientation() eq "C" ) {
        $tbegin = $resultRef->getSubjStart() + $min3;
        $tend = $resultRef->getSubjEnd() - $min5;
      }
      $occBegin{$sName}->{$qName} = $tbegin 
          if ( ! exists $occBegin{$sName}->{$qName} || $occBegin{$sName}->{$qName} > $tbegin );
      $occEnd{$sName}->{$qName} = $tend 
          if ( ! exists $occEnd{$sName}->{$qName} || $occEnd{$sName}->{$qName} < $tend );
    }

    my $cntClipped = 0;
    foreach my $consID ( keys(%{$consRecs}) ) {
      next if ( $consRecs->{$consID}->{'buffer'} == 1 );
      my %alignednr = ();
      foreach my $id (keys %{$occBegin{$consID}}) {
        for (my $i = $occBegin{$consID}->{$id}; $i <= $occEnd{$consID}->{$id}; ++$i) {
          ++$alignednr{$i};
        }
      }
      my $zcount5 = $consRecs->{$consID}->{'leftHPad'};
      my $zcount3 = $consRecs->{$consID}->{'rightHPad'};
      my $conslen = length($consRecs->{$consID}->{'seq'});
      my ($clip5,$clip3,$clipq5,$clipq3) = (0,0,0,0);
      if ($prunecutoff) {
        my $nr5 = $zcount5 + 1;
        while (!$alignednr{$nr5} ||
               $alignednr{$nr5} <= $prunecutoff) { #first base after the Zs
          ++$nr5;
          ++$clip5;
          if ($conslen - $clip5 < 25 ) {
            # this message needs to go to STDOUT (rather than die's STDERR) to be captured by things like $info = `alignAndCallConsensus.pl`;
            print "No or too little alignment was left after pruning fewer than $prunecutoff aligned sequences from the beginning.\n";
            die "Stopped after pruning.\n";
          }
        }
        my $nr3 = $conslen - $zcount3;
        while (!$alignednr{$nr3} ||
               $alignednr{$nr3} <= $prunecutoff) {
          --$nr3;
          ++$clip3;
          if ($conslen - $clip3 < 25 ) {
            print "No or too little alignment was left after pruning fewer than $prunecutoff aligned sequences from the end.\n";
            die "Stopped after pruning.\n";
          }
        }
      }
      if ( $clip5 || $clip3 ) {
        my $tSeq = $consRecs->{$consID}->{'seq'};
        if ( length($tSeq) - $clip5 - $clip3 < 25 ) {
        }
        $tSeq =~ s/^\w{$clip5}// if $clip5;
        $tSeq =~ s/\w{$clip3}$// if $clip3;
        $consRecs->{$consID}->{'seq'} = $tSeq;
        $cntClipped++;
      }
    }
    if ( $cntClipped ) {
      while ( -e "$conFile.$backupIdx" ) {
        $backupIdx++;
      }
      rename "$conFile", "$conFile.$backupIdx";
      open OUT, ">$conFile" or die "Could not open up $conFile for writing!\n";
      foreach my $consID ( keys(%{$consRecs}) ) {
        print OUT ">$consID";
        if ( $consRecs->{$consID}->{'class'} ne "" ) {
          print OUT "#" . $consRecs->{$consID}->{'class'};
        }
        if ( $consRecs->{$consID}->{'desc'} ne "" ) { 
          print OUT "  " . $consRecs->{$consID}->{'desc'};
        }
        print OUT "\n" . $consRecs->{$consID}->{'seq'} . "\n"; 
      }
      close OUT;
    }

    # re-search new consensi vs repseq
    unless ( $engine eq "crossmatch" ) {
      # Always make a database as the sequence may have changed
      system(   "$RMBLAST_DIR/makeblastdb -blastdb_version 4 -out $conFile "
            . "-parse_seqids -dbtype nucl -in $conFile >& /dev/null ");
    }
 
    ( $status, $resultCollection ) = $searchEngineN->search();
    if ( $status )
    {
      print STDERR "Search Parameters: " . $searchEngineN->getParameters() . "\n";
      die "\nERROR from search engine (", $? >> 8, ") \n";
    } 
  
    # Sort the results collection
    $resultCollection->sort(
      sub ($$) {
        $_[ 1 ]->getScore() <=> $_[ 0 ]->getScore();
       }
    );
  } # END PRUNNING
 
  #
  # Loop over all consensi
  #
  my $changedCnt = 0;
  my %filterIndices = ();
  foreach my $consID ( sort { $consRecs->{$a}->{'order'} <=> $consRecs->{$b}->{'order'} } keys(%{$consRecs}) ) {
    next if ( $consRecs->{$consID}->{'stable'} == 1 );
    next if ( $consRecs->{$consID}->{'buffer'} == 1 );
    unless ( $options{'quiet'} ) {
      print "-----\n";
      print "Working on consensus: $consID\n";
    }
    my $removedMinDiv = 10000;
    my $removedMaxDiv = 0;
    my $removedDivCount = 0;
    my $removedDupCount = 0;
    my %queryIDUsed = ();
 

    my $finalAlignCnt = 0;
    # Create a new collection for the filtered set so that we may later sort it
    my $newResultCollection = SearchResultCollection->new();
    for ( my $k = 0 ; $k < $resultCollection->size() ; $k++ ) {
      my $resultRef = $resultCollection->get( $k );
      my $sName = $resultRef->getSubjName();
      $sName = $1 if ( $sName =~ /(\S+)\#.*/ );
      if ( $sName eq $consID ) {
        
        if ( $resultRef->getPctDiverge() <= $maxdiv ) {
          # 'onlyBestAlignment' is deprecated
          if ( $options{'filter_overlapping'} || $options{'onlyBestAlignment'} )
          {
            if ( exists $queryIDUsed{$resultRef->getQueryName()} ) {
              my $overlap = 0;
              foreach my $existingRange ( @{$queryIDUsed{$resultRef->getQueryName()}} ) {
                my $eQStart = $existingRange->[0];
                my $eQEnd = $existingRange->[1];
                my $eSStart = $existingRange->[2];
                my $eSEnd = $existingRange->[3];
                my $cQStart = $resultRef->getQueryStart();
                my $cQEnd = $resultRef->getQueryEnd();
                my $cSStart = $resultRef->getSubjStart();
                my $cSEnd = $resultRef->getSubjEnd();
                if ( ($cQStart >= $eQStart && $cQStart <= $eQEnd) ||
                     ($cQEnd >= $eQStart && $cQEnd <= $eQEnd) ||
                     ($cQStart < $eQStart && $cQEnd > $eQEnd) ||
                     ($cSStart >= $eSStart && $cSStart <= $eSEnd) ||
                     ($cSEnd >= $eSStart && $cSEnd <= $eSEnd) ||
                     ($cSStart < $eSStart && $cSEnd > $eSEnd) ) {
                  $overlap = 1;
                  $removedDupCount++;
                  last;
                }
              }
              $filterIndices{$k} = 1;
              next if ( $overlap );
            }else {
              $queryIDUsed{$resultRef->getQueryName()} = [];
            }
            push @{$queryIDUsed{$resultRef->getQueryName()}}, 
                    [$resultRef->getQueryStart(), $resultRef->getQueryEnd(), $resultRef->getSubjStart(), $resultRef->getSubjEnd()];
            $newResultCollection->add($resultRef);
          }else {
            $newResultCollection->add($resultRef);
          }
          $finalAlignCnt++;
        }else {
          $removedDivCount++;
          $removedMinDiv = $resultRef->getPctDiverge() if ( $resultRef->getPctDiverge() < $removedMinDiv);
          $removedMaxDiv = $resultRef->getPctDiverge() if ( $resultRef->getPctDiverge() > $removedMaxDiv);
        }
      } # if ( $resultRef->getSubjName() eq $consID )...
    } # for ( my $k = 0; $k < $resultCollection->size....

    # Sort the remaining results collection and save in sorted
    # order to $consID.out
    $newResultCollection->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName() <=> $_[ 1 ]->getQueryName() ||
        $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart() ||
        $_[ 1 ]->getQueryEnd() <=> $_[ 0 ]->getQueryEnd();
       }
    );
    open OUT,">$outdir/$consID.out" or die "could not open \'$consID.out\' for writing\n";
    for ( my $k = 0 ; $k < $newResultCollection->size() ; $k++ ) {
      my $resultRef = $newResultCollection->get( $k );
      print OUT "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
    }
    close OUT;
    $newResultCollection = undef;
    

    if ( $removedDivCount > 0 && ! $options{'quiet'} ){
      print "# Removed $removedDivCount high divegence alignments ranging from: $removedMinDiv - $removedMaxDiv divergence\n";
    }
    if ( $removedDupCount > 0 && ! $options{'quiet'} ){
      print "# Removed $removedDupCount alignments to that match the same repseq instance\n";
    }

    if ( $finalAlignCnt < 2 ) {
      print "Only $finalAlignCnt repseq members aligned to $consID in this round.  Cannot generate\n" .
            "a new consensus with this many instances\n"; 
      $consRecs->{$consID}->{'stable'} = 1;
      next;
    }
  
    # Creating new ali file
    my $includeRef = "";
    $includeRef = "-i" if ( $options{'include_reference'} );
    system("$FindBin::RealBin/Linup $includeRef -malignOut $outdir/$consID.malign $outdir/$consID.out > $outdir/$consID.ali");
  
    if ( $options{'html'} ){
      system "$FindBin::RealBin/viewMSA.pl -order start -malign $outdir/$consID.malign";
      system "mv MultipleAlignment.html $consID.html";
    }
    unlink "$outdir/$consID.malign" if ( -e "$outdir/$consID.malign" );
  
    if ( $options{'stockholm'} ) {
      system "$FindBin::RealBin/Linup $includeRef $outdir/$consID.out ~/Matrices/linupmatrix -stockholm > $outdir/$consID.stk";
    }
  
    unless ( $options{'quiet'} ) {
      system "$FindBin::RealBin/scoretotal.pl $outdir/$consID.out";
      my $kimura = `$FindBin::RealBin/CntSubst -cg $outdir/$consID.out | grep kimura`;
      $kimura = $1 if ( $kimura =~ /average kimura-subst: (\S+)/ );
      print "Kimura Divergence: $kimura\n";
    }

    #
    #
    #
    my ($newcons, $hAlignLeft, $hAlignRight, $diffStr) = processAliFile("$outdir/$consID.ali", $ext5done, $ext3done);
    
    my $leftHPad = "H"x($consRecs->{$consID}->{'leftHPad'});
    my $rightHPad = "H"x($consRecs->{$consID}->{'rightHPad'});
    if ($newcons eq $consRecs->{$consID}->{'seq'}) {
      print "Changes: **unchanged**\n" unless ( $options{'quiet'} );
      $consRecs->{$consID}->{'stable'} = 1;

      unless ( $options{'quiet'} ) {
        #print ">$consID\n$newcons\n";
        print "-----\n";
      }
      # Add back the H pads for anchoring purposes only
      $newcons = $leftHPad."\n".$newcons."\n".$rightHPad;

    }else {
      print "Changes:\n$diffStr\n" unless ( $options{'quiet'} );
      $changedCnt++;
      $consRecs->{$consID}->{'iteration'}++;

      if ( $options{'interactive'} ) {
        
        print STDERR "s(kip),c(hangeinbetweenHs),x(pandandchange), b(eginexpand) or 5(\'),e(ndexpand) or 3(\'),##-## (range)\n";
        print STDERR "A range only works if the new and old consensus have the same positions at the start and end of the range.\n";
        my $answer = <STDIN>;
        while ($answer !~ /^[scxbelr35q]$/ && $answer !~ /^\d+[-\s]+\d+/ ) {
          print STDERR "Could not process $answer\nType 's','c','e','b','5','e','3' or '#-###' range.\n" .
                       "Typing \"q\" stops the script";
          $answer = <STDIN>;
        }
        die "Quitting. The consensus file $conFile has not been changed.\n" if $answer eq 'q';;
        chomp $answer;
        if ($answer eq 's') {
           # Do not keep changes
           $newcons = $consRecs->{$consID}->{'seq'};
        } else {
 
          if ($answer eq 'c') {
            # Only allow core changes in new consensus ( i.e. do not allow H pad additions in ).
            $newcons =~ s/^\w{$hAlignLeft}(\w+)\w{$hAlignRight}/$leftHPad\n$1\n$rightHPad/ || 
                print STDERR "Error: processing 'c' hAlignLeft=$hAlignLeft hAlignRight=$hAlignRight,\n" .
                             "       leftHPad=$leftHPad rightHPad=$rightHPad newcons=$newcons\n";
          } elsif ( $answer eq 'x' ) {
            # Trim off edge N's and repad
            $newcons =~ s/^N*(\w+)N*/$leftHPad\n$1\n$rightHPad/;
          } elsif ( $answer =~ /^[b5]$/ ) {
            # Trim off 5' N's and any 3' extension and repad
            $newcons =~ s/^N*(\w+)\w{$hAlignRight}/$leftHPad\n$1\n$rightHPad/  || 
                print STDERR "Error: processing 'b5' hAlignRight=$hAlignRight,\n" .
                             "       leftHPad=$leftHPad rightHPad=$rightHPad newcons=$newcons\n";
          } elsif ( $answer =~ /^[e3]$/ ) {
            # Trim off 3' N's and any 5' extension and repad
            $newcons =~ s/^\w{$hAlignLeft}(\w+)N*/$leftHPad\n$1\n$rightHPad/;
          } elsif ( $answer =~ /^(\d+)\-(\d+)/ ) {
            my $begin = $1 - 1 unless $1 == 0; # half-opened
            my $end = $2;
            my $len = $end - $begin;
            $len = (length $newcons) - $begin if $len > (length $newcons);
            my $newfrag = substr($newcons,$begin,$len);
            my $old5 = "";
            if ($begin) {
              $old5 =  substr($consRecs->{$consID}->{'seq'},0,$begin);
            }
            my $old3 = substr($consRecs->{$consID}->{'seq'},$end);
            $newcons = "$old5"."$newfrag"."$old3";
            $newcons =~ s/^H+/$leftHPad\n/;
            $newcons =~ s/H+$/\n$rightHPad/;
          }
        }
      }else {
        $newcons = $leftHPad."\n".$newcons."\n".$rightHPad;
      }
  
      unless ( $options{'quiet'} ) {
        #print ">$consID\n$newcons\n";
        print "-----\n";
      }

      $consRecs->{$consID}->{'seq'} = $newcons;
    }
  
  } # foreach my $consID...

  # If we have multiple consensi we have to generate an additional file
  # containing all the alignments in one file.
  if ( $numRefineableCons > 1 ) { 
    my $newResultCollection = SearchResultCollection->new();
    for ( my $k = 0 ; $k < $resultCollection->size() ; $k++ ) {
      next if ( exists $filterIndices{$k} );
      my $resultRef = $resultCollection->get( $k );
      $newResultCollection->add($resultRef);
    }
    %filterIndices = ();

    $newResultCollection->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName() <=> $_[ 1 ]->getQueryName() ||
        $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart() ||
        $_[ 1 ]->getQueryEnd() <=> $_[ 0 ]->getQueryEnd();
       }
    );
    open OUT,">$outdir/$conFile.out" or die "could not open \'$conFile.out\' for writing\n";
    for ( my $k = 0 ; $k < $newResultCollection->size() ; $k++ ) {
      my $resultRef = $newResultCollection->get( $k );
      print OUT "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
    }
    close OUT;
    $newResultCollection = undef;
  }
 
  $iterations++;

  if ( $changedCnt ) {
    while ( -e "$conFile.$backupIdx" ) {
      $backupIdx++;
    }
    rename "$conFile", "$conFile.$backupIdx";
    open OUT, ">$conFile" or die "Could not open up $conFile for writing!\n";
    foreach my $consID ( keys(%{$consRecs}) ) {
      print OUT ">$consID";
      if ( $consRecs->{$consID}->{'class'} ne "" ) {
        print OUT "#" . $consRecs->{$consID}->{'class'};
      }
      if ( $consRecs->{$consID}->{'desc'} ne "" ) { 
        print OUT "  " . $consRecs->{$consID}->{'desc'};
      }
      print OUT "\n" . $consRecs->{$consID}->{'seq'} . "\n"; 
    }
    close OUT;
  }

  last if ( ! $options{'refine'} || $changedCnt == 0 );

  if ( $iterations >= $maxRefineIterations ) {
    print "Consensus still changing after $maxRefineIterations iterations...giving up.\n" unless ( $options{'quiet'} );
    last;
  }
}
  
# cleanup
foreach my $ext ( "nog", "nsg", "nsi", "nhr", "nin", "nsq", "nsd" ){
  unlink "$outdir/$conFile.$ext" if ( -s "$outdir/$conFile.$ext" );
}



#
# Process the "ali" file produced by Linup and pull out a consensus
# trimmed based on the extension flags provided by the user and the
# 'H' padding found in the reference.  Also pull out the sections of
# the alignment that show some changes, again disregarding changes
# associated with the H-pads no longer being extended.
#
sub processAliFile {
  my $aliFile = shift;
  my $is5ExtDone = shift;
  my $is3ExtDone = shift;

  open ALI,"<$aliFile" or die "\n\nCould not open ali file \'$aliFile\' for reading!\n\n";

  my $newcons = "";
  my $Hleft = 0;
  my $Hright = 0;
  my @consbit = ();
  my $consline = "";
  my $consbit = "";
  my $diffStr = "";
  while (<ALI>) {
    # consensus     1 TCTTCTGATT----GGTTGGTGGTGAGGTAA-------CA-G...
    if (/^consensus\s+\d+\s+(\S+)/) {
      $newcons .= $1; 
      $consbit = $1;  # Save for comparison of sequences
      $newcons =~ tr/-//d;
      $consline = $_;
      chomp $consline;
    # ref:rep       1 HHHHHTGATT----GGTTGGTGGTGAGGTAA-------CA-G...
    } elsif (/^(ref\:\S+\s+\d+\s+)(\S+)/) {
      my $spacelen = length $1;
      my $refbit = $2;
      if ( $refbit =~ /^(H+)[ACGTNRYMKSW-]/ ) {
        $Hleft = length($1);
      }
      if ( $refbit =~ /[ACGTNRYMKSW-](H+)\s*$/ ) {
        $Hright = length($1);
      }
      if ($refbit ne $consbit) {
        if ( $is5ExtDone && /^ref\:\S+\s+(\d+)\s+(H+)[ACGTN]/ ) {
          # We were done with 5'extension previously.  Any alignment
          # to H's in this round are purely for anchoring purposes and
          # should not contribute the consensus.
          $Hleft = length $2;
          $newcons =~ s/^\w{$Hleft}//;
          # Now that we are done with extension only record differences
          # if they occur after the H padding.
          my ($tempcons,$tempref) = ($consbit,$refbit);
          $tempcons =~ s/^\w{$Hleft}//;
          $tempref =~ s/^H+//;
          next if $tempcons eq $tempref;
        } elsif ($is3ExtDone) {
          if ($Hright && /^ref\:\S+\s+\d+\s+(H+)\s/ ) {
            $Hright += length $1;
            next;
          } elsif ( /^ref\:\S+\s+\d+\s+\S*?(H+)\s/ ) {
            $Hright = length $1;
            my ($tempcons,$tempref) = ($consbit,$refbit);
            $tempcons =~ s/\w{$Hright}$//;
            $tempref =~ s/H+$//;
            next if $tempcons eq $tempref;
          }
        }
  
        my $midline = " " x $spacelen;
        my @refbases = split '',$refbit;
        my @consbases = split '',$consbit;
        for (my $i = 0; $i <= $#refbases; ++$i) {
          if ($refbases[$i] eq $consbases[$i]) {
            $midline .= " ";
          } elsif ($refbases[$i] =~ /[CTY]/ && $consbases[$i] =~ /[CTY]/ ||
                   $refbases[$i] =~ /[AGR]/ && $consbases[$i] =~ /[AGR]/ ) {
            $midline .= "i";
          } elsif ($refbases[$i] =~ /[CTY]/ && $consbases[$i] =~ /[AGR]/ ||
                   $refbases[$i] =~ /[AGR]/ && $consbases[$i] =~ /[CTY]/ ) {
            $midline .= "v";
          } elsif ($refbases[$i] =~ /\w/ && $consbases[$i] =~ /-/ ||
                   $refbases[$i] =~ /-/ && $consbases[$i] =~ /\w/ ) {
            $midline .= "-";
          } else {
            $midline .= '?';
          }
        }
        $diffStr .= "$consline\n$midline\n$_\n";
      }
    }
  }
  close ALI;
  $newcons =~ s/(\w){$Hright}$// if $Hright && $is3ExtDone;
 
  return( $newcons, $Hleft, $Hright, $diffStr );
}


#
# Read in the consensus file and build up a data structure for use
# by the main code.
#
sub processConFile {
  my $conFile = shift;

  my($filename, $dirs, $suffix) = fileparse($conFile);
  open CON,"<$conFile" or die "\n\nCould not open consensus file \'$conFile\' for reading!\n\n";
  my $numBuffers = 0;
  my $idx = 1;
  my %consRecs = ();
  my $id;
  my $seq;
  my $desc;
  my $class;
  while (<CON>) {
    if ( /^>(\S+)\s+(.*)$/ ) {
      my $tID = $1;
      my $tDesc = $2;
      $class = "";
      if ( $seq ) {
        my $isBuffer = 0;
        if ( $id =~ /(\S+)\#(\S+)/ ) {
          $id = $1;
          $class = $2;
          if ( $class =~ /buffer/i ) {
            $numBuffers++;
            $isBuffer = 1;
          }
        }
        if ( exists $consRecs{$id} ) {
          die "\n\nConsensus file contains a duplicate ID = \'$id\'.\n\n";
        }
        if ( $id eq $filename ) {
          die "\n\nConsensus file contains a identifier with the same name as the file $id!\n\n";
        }
        # LEGACY "Z" padding.
        if ( $seq =~ /^Z+/ || $seq =~ /.*Z+$/ ) {
          $seq =~ s/Z/H/g;
        }
        my $leftHPad = 0;
        my $rightHPad = 0;
        if ( $seq =~ /^(H+)(\S+?)(H*)$/ ) {
          $leftHPad = length($1);
          $rightHPad = length($3) if ( $3 );
          $seq = $2;
        }elsif ( $seq =~ /^(\S+?)(H*)$/ ) {
          $rightHPad = length($2) if ( $2 );
          $seq = $1;
        }
        if ( $seq !~ /^(H*)[ACGTRYMKSWN]+(H*)$/ )
        {
          my $problems = "";
          while ( $seq =~ s/(([ACGTRYMKSWN]{0,3})([^ACGTRYMKSWN]+)([ACGTRYMKSWN]{0,3}))//g ) {
            $problems .= "\'" . lc($2) . "$3". lc($4) . "\', ";
          }
          die "\n\nConsensus $id contains characters that are not supported.  Currently\n" .
              "this program supports the use of ACGT an the IUB codes RYMKSWN.  Please\n" .
              "correct the sequence and rerun. Problems: $problems\n\n";
        }
        $consRecs{$id} = { 'seq' => $seq, 
                           'class' => $class,
                           'desc' => $desc,
                           'buffer' => $isBuffer,
                           'iteration' => 0,
                           'stable' => 0,
                           'leftHPad' => $leftHPad,
                           'rightHPad' => $rightHPad,
                           'order' => $idx };
        $idx++;
      }
      $id = $tID;
      $desc = $tDesc;
      $seq = "";
      next;
    }
    s/[\n\r\s]//g;
    $seq .= uc($_);
  }
  close CON;
  
  # Trailing case
  if ( $seq ) {
    my $isBuffer = 0;
    if ( $id =~ /(\S+)\#(\S+)/ ) {
      $id = $1;
      $class = $2;
      if ( $class =~ /buffer/i ) {
        $numBuffers++;
        $isBuffer = 1;
      }
    }
    if ( exists $consRecs{$id} ) {
      die "\n\nConsensus file contains a duplicate ID = \'$id\'.\n\n";
    }
    if ( $id eq $filename ) {
      die "\n\nConsensus file contains a identifier with the same name as the file $id!\n\n";
    }
    my $leftHPad = 0;
    my $rightHPad = 0;
    if ( $seq =~ /^(H+)(\S+?)(H*)$/ ) {
      $leftHPad = length($1);
      $rightHPad = length($3) if ( $3 );
      $seq = $2;
    }elsif ( $seq =~ /^(\S+?)(H*)$/ ) {
      $rightHPad = length($2) if ( $2 );
      $seq = $1;
    }
    if ( $seq !~ /^(H*)[ACGTRYMKSWN]+(H*)$/ )
    {
      my $problems = "";
      while ( $seq =~ s/(([ACGTRYMKSWN]{0,3})([^ACGTRYMKSWN]+)([ACGTRYMKSWN]{0,3}))//g ) {
        $problems .= "\'" . lc($2) . "$3". lc($4) . "\', ";
      }
      die "\n\nConsensus $id contains characters that are not supported.  Currently\n" .
          "this program supports the use of ACGT an the IUB codes RYMKSWN.  Please\n" .
          "correct the sequence and rerun. Problems: $problems\n\n";
    }
    $consRecs{$id} = { 'seq' => $seq, 
                       'class' => $class,
                       'desc' => $desc,
                       'buffer' => $isBuffer,
                       'iteration' => 0,
                       'stable' => 0,
                       'leftHPad' => $leftHPad,
                       'rightHPad' => $rightHPad,
                       'order' => $idx };
    $idx++;
  }
 
  my $numRefineableCons = scalar(keys(%consRecs))-$numBuffers;

  return( $numRefineableCons, $numBuffers, \%consRecs );
}

1;
