#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) alignAndCallConsensus.pl
##  Author:
##      Arian Smit         asmit@systemsbiology.org
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A script to align a consensus to a set of putative
##      family instances and recall the consensus based on the
##      linup of all alignments.
##

=head1 NAME

alignAndCallConsensus.pl - Align TE consensus to a set of instances and call consensus

=head1 SYNOPSIS

  alignAndCallConsensus.pl [-c(onsensus) <*.fa>] [-e(lements) <*.fa>] 
                 [-d(ivergencemax) #] [-sc(ore) #] [-mi(nmatch) #]
                 [-h(tml)] [-st(ockholm)] [-q(uiet)]
                 [-b(andwidth) #] [-cr(ossmatch)] [-rm(blast)]
                 [-re(fine)|-re(fine) #> [-fi(lter_overlap)]
                 [-f(inishedext) <str>] [-ma(trix) <str>]
                 [-inc(lude_reference) <val>]
                 [-p(rune) <val>] [-qu(oteprune) <val>]
                 [-int(eractive)] [-help]

  Example:

    ./alignAndCallConsensus.pl

        Assumes that a consensus file named "rep" and an elements file 
        named "repseq" ( representative sequences ) exists in the 
        current working directory.  All defaults will be used to
        generate a single transitive multiple alignment and a new 
        consensus.

  or

    ./alignAndCallConsensus.pl -c mycons.fa -e myinstances.fa


  This tool was original written by Arian Smit as two different scripts known 
  affectionately as "dothemsimple.pl" and "dothemmultiple.pl".  As part of a 
  large set of small command line tools this script performs an alignment of 
  a putative TE family consensus or a set of related TE family consensi to a 
  representative set of instances for the family(ies).  The alignments are then
  used to generate a transitively induced multiple alignment ( aka "line up")
  for the purposes of consensus calling.  The new consensus may then be used 
  to realign to the instances and another consensus may be called.  In a sense 
  this is a hill climbing optimisation.  The output consensus may be better 
  than the input consensus but it still may not be optimal or even stable.  By 
  applying this method iteratively the consensus will improve and stabilize.

  The tool was written to support multiple use cases:

  Consensus Extension:
     Provided a consensus from a de-novo TE discovery tool, this tool may
     be used to manually extend a consensus sequence.  The "H" IUB code
     has been repurposed for this task.  By prefixing/suffxing "H" characters
     to the input consensus, flanking bases may be pulled into the final 
     alignment. The matrices used by this utility score any match to an
     "H" positively (+3).  This has the effect of extending the alignment
     by the number of "H" on either end of the alignment.  This assumes that
     the -elements file contains flanking sequence in advance.
     NOTE: This will repad the final consensus with the same number of H's
           that it started with to assist with continuing the extension.

       E.g.

          initial consensus:    HHHHAACGTTAGCGGATTACGAGGGCAHHHH
          elements:             AATAAACGTTAGCGGATTACGAGGGCAGATA
                                ACTAA-CGTAAGCGAATTACGA-GGCAGCTC
                                ACTAA-CGTTAGCGAATTACGAGGGCAGTTG
                                ACTAA-CGTTAGCGAATTACGAGGGCAGGTT

          new consensus:        ACTAA-CGTTAGCGAATTACGAGGGCAGNTN
          updated positions:    **** *        *            * *  


  Iterative Refinement:
     The initial consensus may be rudimentary or simply derived from a
     different representative set.  In this case the new consensus
     output at the end of the the run may differ from the original.  In
     this case this tool may be re-run on the updated consensus in an
     iterative fashion until the input consensus matches the output
     conensus.


  Subfamily Consensus Building:
     After performing a subfamily analysis using tools such as Coseg,
     one often wants to build a formal alignment and consensus for each
     of a set of related subfamilies.  Typically all starting consensi
     a used to gather a combined set of instances.  The consenus
     file contains all of the starting subfamilies.  This script will
     align each instance to all consensi and assign each one to it's
     best match. Then each alignment set will be used to calculate a new 
     cosensus for each family in the same fashion as the single consensus
     case.


  Output:
     The program saves the original input consensus file as <filename>.#, where 
     the number begins with '1' and is incremented up until 25 at which time
     the program dies and requires backups be removed before proceeding. The new
     consensus is saved as the original input consensus file name.


=head1 DESCRIPTION

The options are:

=over 4

=item -c(onsensus) <*.fa> 

A FASTA file containing a single consensus sequence for a TE family. In 
addition to the standard nucleotide characters the consensus may contain
a limited IUB code set: N, R, Y, M, K, S, and W.  The weakly specified 
IUB codes W,H,B,V and D are reserved for other uses.  

The consensus may be optionally padded on either side (5'/3') with the 
special purpose "H" pad character.  Using a special matrices which score
positively any match to an "H", this utility is able to pull in flanking
sequence into the alignments.  Default: If this option is not provided
the file "rep" in the current working directory is assumed to be the
consensus file.

=item -e(lements) <*.fa>

A multi-sequence FASTA file containing sequences believed to be elements of 
a single TE family.  Sequences may (and often should) contain flanking sequence
for possible extension using the "H" padding character (see consensus parameter).
Default: If this option is not provided the file "repseq" in the current
working directory is assumed to be the elements file.


=item -d(ivergencemax) #
 
Filters aligned element sequences above this % divergence. Default: 60 

=item -sc(ore) #

The minimum aligned score. Default: 200

=item -mi(nmatch) #

The minimum word size for the alignment algorithm. Default: 7

=item -b(andwidth) #

The alignment bandwidth for cross_match ( or equivalent x-drop cutoffs for rmblast )
Default: 40

=item -ma(trix) <value>

Default 25; If a number 14,18,20 or 25 is given the matrix ##p41g.matrix is chosen.
If given 14p35 (37,39,41,43,45,47,49,51,53) #####g.matrix is chosen
 
=item -f(inishedext) <value>

Default off; given N 'H's on an end of the reference, the program wil extend 
the consensus by N bases and pad the new consensus with N 'H's.
-f 5 : 5' extension done (the 5' end won't grow nor get 'H's attached)
-f 3 : 3' extension done
-f b : both extensions done

=item -inc(lude_reference)

Include the supplied consensus sequence in the multiple alignment.  Used in cases
were the starting consnesus is a single instance of the family.

=item -p(rune) # or "#n#"

Prunes the consensus sequence on either side to the point that more than # 
sequences are aligned. If the string form "#n#" is used ....

=item -qu(oteprune) #

Same as prune except that instead of pruning to the point that more than #
sequences are aligned it prunes until the number of aligned sequences increases
by >= # times.  Only checks the first and last 100bp and when 3 or more seqs
join.

=item -h(tml)       

Generate an HTML visualization of the alignment using viewMSA.pl. 

=item -st(ockholm)     

Creates a Stockholm format (.stk) file on top of the normal Linup alignment file

=item -rm(blast)

Use RMBlast instead of cross_match

=item -re(fine) or -re(fine) <value>

Iterate the process of aligning the derived consensus to the set
of instance sequences.  This is a simple optimisation step which 
(like EM) which typically stabilizes after a small number of
iterations.  The default is to attempt 5 iterations but a higher
value may be specified with this option.

=item -fi(lter_overlap)

When a single sequence in the elements file produces multiple sub
alignments to the consensus, this option will remove any sub alignment
which overlaps the same consensus range *or* the same sequence
range as a higher scoring sub alignment.

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
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
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

usage() if ( $options{'help'} );

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
    print "\n\nConsensus file \"$conFile\" is missing or empty! Use '-help' for option help.\n\n";
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

my $elesFile = "repseq";
if ( exists $options{'elements'} ) {
  if ( -s $options{'elements'} ) {
    $elesFile = $options{'elements'}; 
  }else {
    print "\n\nElements file \"$elesFile\" is missing or empty!\n\n";
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
$minscore = $options{'minscore'} if ( exists $options{'minscore'} );

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
        print OUT ">$consID\n" . $consRecs->{$consID}->{'seq'} . "\n"; 
      }
      close OUT;
    }

    # re-search new consensi vs repseq
    unless ( $engine eq "crossmatch" ) {
      # Always make a database as the sequence may have changed
      system(   "$RMBLAST_DIR/makeblastdb -blastdb_version 4 -out $conFile "
            . "-parse_seqids -dbtype nucl -in $conFile >& /dev/null ");
    }
 
    #print STDERR "Search Parameters: " . $searchEngineN->getParameters() . "\n";
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
  foreach my $consID ( keys(%{$consRecs}) ) {
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
 
    open OUT,">$outdir/$consID.out" or die "could not open \'$consID.out\' for writing\n";

    my $finalAlignCnt = 0;
    for ( my $k = 0 ; $k < $resultCollection->size() ; $k++ ) {
      my $resultRef = $resultCollection->get( $k );
      if ( $resultRef->getSubjName() eq $consID ) {
        
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
              next if ( $overlap );
            }else {
              $queryIDUsed{$resultRef->getQueryName()} = [];
            }
            push @{$queryIDUsed{$resultRef->getQueryName()}}, 
                    [$resultRef->getQueryStart(), $resultRef->getQueryEnd(), $resultRef->getSubjStart(), $resultRef->getSubjEnd()];
            print OUT "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
          }else {
            print OUT "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
          }
          $finalAlignCnt++;
        }else {
          $removedDivCount++;
          $removedMinDiv = $resultRef->getPctDiverge() if ( $resultRef->getPctDiverge() < $removedMinDiv);
          $removedMaxDiv = $resultRef->getPctDiverge() if ( $resultRef->getPctDiverge() > $removedMaxDiv);
        }
      } # if ( $resultRef->getSubjName() eq $consID )...
    } # for ( my $k = 0; $k < $resultCollection->size....
    if ( $removedDivCount > 0 && ! $options{'quiet'} ){
      print "# Removed $removedDivCount high divegence alignments ranging from: $removedMinDiv - $removedMaxDiv divergence\n";
    }
    if ( $removedDupCount > 0 && ! $options{'quiet'} ){
      print "# Removed $removedDupCount alignments to that match the same repseq instance\n";
    }
    close OUT;

    if ( $finalAlignCnt < 2 ) {
      print "Only $finalAlignCnt repseq members aligned to $consID in this round.  Cannot generate\n" .
            "a new consensus with this many instances\n"; 
      $consRecs->{$consID}->{'stable'} = 1;
      next;
    }
  
    # Creating new ali file
    my $includeRef = "";
    $includeRef = "-i" if ( $options{'include_reference'} );
    system("$FindBin::RealBin/Linup $includeRef $outdir/$consID.out > $outdir/$consID.ali");
  
    if ( $options{'html'} ){
      system "$FindBin::RealBin/viewMSA.pl -order start -malign $outdir/out.malign";
      system "mv MultipleAlignment.html $consID.html";
    }
  
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
    
    if ($newcons eq $consRecs->{$consID}->{'seq'}) {
      print "Changes: **unchanged**\n" unless ( $options{'quiet'} );
      $consRecs->{$consID}->{'stable'} = 1;

      unless ( $options{'quiet'} ) {
        print ">$consID\n$newcons\n";
        print "-----\n";
      }

    }else {
      print "Changes:\n$diffStr\n" unless ( $options{'quiet'} );
      my $leftHPad = "H"x($consRecs->{$consID}->{'leftHPad'});
      my $rightHPad = "H"x($consRecs->{$consID}->{'rightHPad'});
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
        print ">$consID\n$newcons\n";
        print "-----\n";
      }

      $consRecs->{$consID}->{'seq'} = $newcons;
    }
  
  } # foreach my $consID...

  $iterations++;

  if ( $changedCnt ) {
    while ( -e "$conFile.$backupIdx" ) {
      $backupIdx++;
    }
    rename "$conFile", "$conFile.$backupIdx";
    open OUT, ">$conFile" or die "Could not open up $conFile for writing!\n";
    foreach my $consID ( keys(%{$consRecs}) ) {
      print OUT ">$consID\n" . $consRecs->{$consID}->{'seq'} . "\n"; 
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
unlink "$outdir/out.malign" if ( -e "$outdir/out.malign" );



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

  open CON,"<$conFile" or die "\n\nCould not open consensus file \'$conFile\' for reading!\n\n";
  my $numRefineableCons = 0;
  my $numBuffers = 0;
  my %consRecs = ();
  my $id;
  my $seq;
  while (<CON>) {
    if ( /^>(\S+)/ ) {
      my $tID = $1;
      if ( $seq ) {
        my $isBuffer = 0;
        if ( $id =~ /.*\#buffer/i ) {
          $numBuffers++;
          $isBuffer = 1;
        }else {
          $numRefineableCons++;
        }
        if ( exists $consRecs{$id} ) {
          die "\n\nConsensus file contains a duplicate ID = \'$id\'.\n\n";
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
                           'buffer' => $isBuffer,
                           'iteration' => 0,
                           'stable' => 0,
                           'leftHPad' => $leftHPad,
                           'rightHPad' => $rightHPad };
      }
      $id = $tID;
      $seq = "";
      next;
    }
    s/[\n\r\s]//g;
    $seq .= uc($_);
  }
  if ( $seq ) {
    my $isBuffer = 0;
    if ( $id =~ /.*\#buffer/i ) {
      $numBuffers++;
      $isBuffer = 1;
    }else {
      $numRefineableCons++;
    }
    if ( exists $consRecs{$id} ) {
      die "\n\nConsensus file contains a duplicate ID = \'$id\'.\n\n";
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
                       'buffer' => $isBuffer,
                       'iteration' => 0,
                       'stable' => 0,
                       'leftHPad' => $leftHPad,
                       'rightHPad' => $rightHPad };
  }
  close CON;

  return( $numRefineableCons, $numBuffers, \%consRecs );
}

1;
