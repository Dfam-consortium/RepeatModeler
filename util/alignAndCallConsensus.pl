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

  alignAndCallConsensus.pl -c(onsensus) <*.fa> -e(lements) <*.fa> | -de(faults)
                 [-d(ivergencemax) #] [-sc(ore) #] [-mi(nmatch) #]
                 [-b(andwidth) #] [-cr(ossmatch)] [-rm(blast)]
                 [-f(inishedext) <str>] [-ma(trix) <str>]
                 [-ht(ml)] [-st(ockholm)] [-q(uiet)]
                 [-re(fine)|-re(fine) #] [-fi(lter_overlap)]
                 [-inc(lude_reference)] [-outdir <val]
                 [-p(rune) <val>] [-int(eractive)] 
                 [-hp(ad) #] [-h(elp)]

  Example:

    ./alignAndCallConsensus.pl -defaults

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
     In some cases the starting consensus may be fairly distant from
     the sequences being aligned.  In this case the alignment and the
     calculated consensus may still be suboptimal after one iteration.
     The script may be called again using the new consensus to iteratively
     improve the alignment and consensus or the -refine option may be used
     automatically reiterate until the consensus stabilizes (doesn't change)
     or a maximum iteration count is reached.

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

  Buffer Sequences:


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

Default off; Indicate that one or both ends are not to be extended even if they 
have 'H's prefixed/suffixed to the sequence.
-f 5 : 5' extension done, do not incorporate bases matching the H-pad into consensus.
-f 3 : 3' extension done, ...
-f b : both extensions done, ...

=item -inc(lude_reference)

Include the supplied consensus sequence in the multiple alignment.  Used in cases
were the starting consnesus is a single instance of the family.

=item -p(rune) # or "#n#"

Prunes the consensus sequence on either side to the point that more than # 
sequences are aligned. If the string form "#n#" is used ....

=item -ht(ml)       

Generate an HTML visualization of the alignment using viewMSA.pl. 

=item -hp(ad) #

Append/prepend 'H' characters to the sequence so that there are the given number on each
side.  The current number of 'H's will be taken into account so the actual number may grow
or shrink the sequence.  This is performed prior to the first search.

=item -st(ockholm)     

Creates a Stockholm format (.stk) file on top of the normal Linup alignment file

=item -rm(blast)

Use RMBlast instead of cross_match

=item -re(fine) or -re(fine) <value>

Iterate the process of aligning the derived consensus to the set
of instance sequences.  This is a simple optimisation step which 
(like EM) typically stabilizes after a small number of
iterations.  The default, if a value is not passed after '-refine' 
is to attempt up to five iterations but a higher/lower value may 
be specified with this option.

=item -fi(lter_overlap)

When a single sequence in the elements file produces multiple sub
alignments to the consensus, this option will remove any sub alignment
which overlaps the same consensus range *or* the same sequence
range as a higher scoring sub alignment.


=back

=head1 SEE ALSO

ReapeatModeler, dothemsimple.pl, extendcons.pl, createxmoutstk.pl etc..

=head1 COPYRIGHT

Copyright 2019-2021 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR(s)

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
my $Version = $RepModelConfig::VERSION;

# A global hash to hold program timers
my %TimeBefore = ();


#
# Paths
#
my $phrapDir = "/usr/local/phrap";
my $matrixDir = "$FindBin::RealBin/../Matrices";
my $RMBLAST_DIR = $RepModelConfig::configuration->{'RMBLAST_DIR'}->{'value'};
my $defaultEngine = "rmblast";

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
                    '-divergencemax|d=s', 
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
                    '-html|ht',
                    '-hpad|hp=i',
                    '-quiet|q',
                    '-refine|re:s',
                    '-bandwidth|b=s',
                    '-crossmatch|cr',
                    '-rmblast|rm',
                    '-filter_overlapping|fi',
                    '-debug',
                    '-onlyBestAlignment', # Deprecated, now 'filter_overlapping'
);
my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

# TODO: NOT IMPLEMENTED
# '-quoteprune|qu=s',
#
# -qu(oteprune) #
#
# Same as prune except that instead of pruning to the point that more than #
# sequences are aligned it prunes until the number of aligned sequences increases
# by >= # times.  Only checks the first and last 100bp and when 3 or more seqs
# join.
#

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0 | less";
  exit( 1 );
}

unless ( ! $options{'help'} && (
         $options{'defaults'} || ( $options{'consensus'} && $options{'elements'} )))
{
   print "\n\nMissing either '-defaults' or '-consensus <file> -elements <file>'\nUse '$0 -h' to view the help.\n\n";
   exit(1);
}

if ( exists $options{'refine'} && exists $options{'interactive'} ) {
   print "\n\nOptions -refine and -interactive are mutually exclusive!\nUse '$0 -h' to view the help.\n\n";
   exit(1);
}

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
my $hpadCnt;
$hpadCnt = $options{'hpad'} if ( $options{'hpad'} );
my ( $numRefineableCons, $numBuffers, $consRecs, $seqsWithHpads ) = &processConFile( $conFile, $hpadCnt );

# If -hpad is provided, it is more than likely that the consensus sequences have
# been modified to include a different number of H's.  We should right out the
# consensus file before the search.
if ( $options{'hpad'} ) {
  open OUT, ">$conFile" or die "Could not open up $conFile for writing!\n";
  foreach my $consID ( keys(%{$consRecs}) ) {
    print OUT ">$consID";
    if ( $consRecs->{$consID}->{'class'} ne "" ) {
      print OUT "#" . $consRecs->{$consID}->{'class'};
    }
    if ( $consRecs->{$consID}->{'desc'} ne "" ) { 
      print OUT "  " . $consRecs->{$consID}->{'desc'};
    }
    print OUT "\n" . "H"x($options{'hpad'}) . $consRecs->{$consID}->{'seq'} . "H"x($options{'hpad'}) . "\n"; 
  }
  close OUT;
}

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
if ( exists $options{'crossmatch'} ){
  if ( $phrapDir eq "" || ! -x "$phrapDir/cross_match"  ) {
    print "\n\nThis option only works if the \$phrapDir variable near the top of the script is\n";
    print "defined.  This is currently a non-supported feature.\n";
    exit(1);
  }
  $engine = "crossmatch"
}
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
elapsedTime('main');
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
  if ( $prunecutoff > 0 ) {
    print "# Prune alignment edges where coverage is less than $prunecutoff sequences.\n";
  }
  if ( $seqsWithHpads ) {
    print "# Extension Mode, $seqsWithHpads sequences have Hpads\n";
  }
  if ( exists $options{'refine'} ) {
    print "# Refinement Mode,  max $maxRefineIterations iterations\n";
  }
  if ( $firstFreeBackupIdx > 1 ) {
    print "# Starting Round Index: $firstFreeBackupIdx\n";
  }
  print "#--------------------------------------------------------------------\n";
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
            . "-parse_seqids -dbtype nucl -in $conFile >/dev/null 2>&1");
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

  print "#\n";
  print "# ITERATION " . ( $iterations + 1 ) . "\n";
  print "#\n";
  print "#--------------------------------------------------------------------\n";

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
      my $subjName = $resultRef->getSubjName();
      my $sName = $resultRef->getSubjName();
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
      # Subject is the reference so these are reference coordinates
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

      #print "Coverage:\n";
      #foreach my $k ( sort{$a <=> $b} keys(%alignednr) ) {
      #  print ", $k:" . $alignednr{$k};
      #}
      #print "\n";

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
        my $nr3 = $conslen + $zcount5;
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
      print "Consensus Pruned\n";
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
            . "-parse_seqids -dbtype nucl -in $conFile >/dev/null 2>&1 ");
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
  my $numSkipped = 0;
  foreach my $consID ( sort { $consRecs->{$a}->{'order'} <=> $consRecs->{$b}->{'order'} } keys(%{$consRecs}) ) {
    next if ( $consRecs->{$consID}->{'stable'} == 1 );
    next if ( $consRecs->{$consID}->{'buffer'} == 1 );
    print "Working on $consID\n";
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
      system "$FindBin::RealBin/Linup $includeRef $outdir/$consID.out -matrix $FindBin::RealBin/../Matrices/linupmatrix -stockholm > $outdir/$consID.stk";
    }
  
    unless ( $options{'quiet'} ) {
      &scoretotal("$outdir/$consID.out");
      #out.CGmodified: average kimura-subst: 0.149304766710626
      my $kimura = `$FindBin::RealBin/CntSubst -cg $outdir/$consID.out | grep kimura`;
      $kimura = $1 if ( $kimura =~ /average kimura-subst: (\S+)/ );
      #out.CGmodified: total aligned length: 15249 bp
      my $aligned_bp = `$FindBin::RealBin/CntSubst -cg $outdir/$consID.out | grep aligned`;
      $aligned_bp = $1 if ( $aligned_bp =~ /length: (\S+) bp/ );
      print "Kimura Divergence: " . sprintf("%0.3f",$kimura*100) . " % ( $aligned_bp aligned bps )\n";
    }

    my $leftHPad = "H"x($consRecs->{$consID}->{'leftHPad'});
    my $rightHPad = "H"x($consRecs->{$consID}->{'rightHPad'});

    my ( $consHLeft, $consCore, $consHRight, $consAligned, $refAligned, $diffStr ) = processAliFile("$outdir/$consID.ali", $ext5done, $ext3done);
    #print "ALI CONS: $consHLeft ... $consCore ... $consHRight\n";

    my $testCons = $consCore;
    $testCons = $consHLeft . $testCons unless( $ext5done );
    $testCons .= $consHRight unless( $ext3done );
    my $newcons = $testCons;

    if ($testCons eq $consRecs->{$consID}->{'seq'}) {
      print "Consensus Changes: **unchanged**\n" unless ( $options{'quiet'} );
      $consRecs->{$consID}->{'stable'} = 1;

      unless ( $options{'quiet'} ) {
        print "#--------------------------------------------------------------------\n";
      }
    }else {
      print "Consensus Changes:\n\n$diffStr\n" unless ( $options{'quiet'} );
      #print "new: $testCons\n";
      #print "old: " . $consRecs->{$consID}->{'seq'} ."\n";
      $changedCnt++;
      $consRecs->{$consID}->{'iteration'}++;

      my $intMessage = "";
      if ( $options{'interactive'} ) {
        
        print STDERR "s(kip),c(hangeinbetweenHs),x(pandandchange), b(eginexpand) or 5(\'),e(ndexpand) or 3(\'),##-## (cons range),d(one)\n";
        my $answer = <STDIN>;
        while ($answer !~ /^[scxbelr35d]$/ && $answer !~ /^\d+[-\s]+\d+/ ) {
          print STDERR "Could not process $answer\nType 's','c','e','b','5','e','3','#-###' range, or 'd'.\n" .
                       "Typing \"d\" quits early keeping only the changes made so far.\n";
          $answer = <STDIN>;
        }
        chomp $answer;
        if ( $answer eq "d" ) {
          print "Done! Consensus file ($conFile) has been updated with any previously made selections.\n";
          &saveNewCons( $conFile, $consRecs );
          &cleanup( $outdir, $conFile );
          exit();
        }
        if ($answer eq 's') {
           # Do not keep changes - take the previous consensus 
           $numSkipped++;
           $newcons = $consRecs->{$consID}->{'seq'};
           if ( $numRefineableCons == $numSkipped ) {
             print "Changes skipped for all consensi. Done!\n";
             &cleanup( $outdir, $conFile );
             exit();
           }
        } else {
          #print "hAlignLeft = $hAlignLeft, hAlignRight=$hAlignRight, leftHPad=$leftHPad, rightHPad=$rightHPad\n";
          if ($answer eq 'c') {
            # Only allow core changes in new consensus ( i.e. do not allow H pad additions in ).
            #print "processing 'c' hAlignLeft=$hAlignLeft hAlignRight=$hAlignRight,\n" .
            #      "       leftHPad=$leftHPad rightHPad=$rightHPad cons=$newcons\n";
            #$newcons =~ s/^\w{$hAlignLeft}(\w+)\w{$hAlignRight}/$1/ || 
            #    print STDERR "Error: processing 'c' hAlignLeft=$hAlignLeft hAlignRight=$hAlignRight,\n" .
            #                 "       leftHPad=$leftHPad rightHPad=$rightHPad newcons=$newcons\n";
            $newcons = $consCore;
            $intMessage = "Keeping only core changes, ignoring the sequence in the H-pad regions.";
          } elsif ( $answer eq 'x' ) {
            # Trim off edge N's and repad
            $newcons = $consHLeft . $consCore . $consHRight;
            $newcons =~ s/^N*(\w+)N*/$1/;
            $intMessage = "Keeping both 5\' and 3\' H-pad changes.";
          } elsif ( $answer =~ /^[b5]$/ ) {
            # Trim off 5' N's and any 3' extension and repad
            $newcons = $consHLeft . $consCore; 
            $newcons =~ s/^N*(\w+)/$1/  || 
                print STDERR "Error: processing 'b5' consHLeft=$consHLeft,\n" .
                             "       leftHPad=$leftHPad rightHPad=$rightHPad newcons=$newcons\n";
            $intMessage = "Keeping only 5\' H-pad changes.";
          } elsif ( $answer =~ /^[e3]$/ ) {
            # Trim off 3' N's and any 5' extension and repad
            $newcons = $consCore . $consHRight; 
            $newcons =~ s/^(\w+)N*/$1/;
            $intMessage = "Keeping only 3\' H-pad changes.";
          } elsif ( $answer =~ /^(\d+)\-(\d+)/ ) {
            $intMessage = "Keeping only consensus range $1-$2.";
            my $rangeStart = $1;
            my $rangeEnd = $2;
            my $rangeColStart = 0;
            my $rangeColEnd = 0;
            my $baseCnt = 0;
            for ( my $i = 0; $i < length($consAligned); $i++ ){
              $baseCnt++ if ( substr($consAligned,$i,1) ne "-" );
              $rangeColStart = $i;
              last if ( $baseCnt == $rangeStart );
            }
            if ( $baseCnt != $rangeStart ) { print "Range Error: could not find start position of $rangeStart!\n"; }
            for ( my $i = $rangeColStart+1; $i < length($consAligned); $i++ ){
              $baseCnt++ if ( substr($consAligned,$i,1) ne "-" );
              $rangeColEnd = $i;
              last if ( $baseCnt == $rangeEnd );
            }
            if ( $baseCnt != $rangeEnd ) { print "Range Error: could not find end position of $rangeStart!\n"; }
            #print "Identified range columns: $rangeColStart - $rangeColEnd\n";
            #print "Cons for this range: " . substr($consAligned,$rangeColStart, ($rangeColEnd - $rangeColStart + 1)) . "\n";
            #print "Ref for this range: " . substr($refAligned,$rangeColStart, ($rangeColEnd - $rangeColStart + 1)) . "\n";
            substr($refAligned,$rangeColStart, ($rangeColEnd - $rangeColStart + 1)) = substr($consAligned,$rangeColStart, ($rangeColEnd - $rangeColStart + 1));
            $newcons = $refAligned;
            $newcons =~ s/-//g;
            $newcons =~ s/^H*(.*)/$1/;
            $newcons =~ s/([^H]*)H+$/$1/;
            #print "newCons = $newcons\n";
          }
        }
        if ( $options{'debug'} ) {
          print "Newcons: $leftHPad$newcons$rightHPad\n";
        }
      } # interactive
  
      if ( $intMessage ne "" ) {
        print "$intMessage\n";
      }

      unless ( $options{'quiet'} ) {
        print "#--------------------------------------------------------------------\n";
      }

      $consRecs->{$consID}->{'seq'} = $newcons
    } # has changes
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
    &saveNewCons($conFile,$consRecs);
  }

  if ( exists $options{'refine'} ) {
    last if ( $changedCnt == 0);
    if ( $iterations >= $maxRefineIterations ) {
      if ( $changedCnt != 0 ) {
        print "#\n";
        print "# WARN: Consensus still changing after $maxRefineIterations iterations. May need to continue refinement.\n" unless ( $options{'quiet'} );
        print "#\n";
      }
      last;
    }
  }elsif ( $options{'interactive'} ) {
    last if ( $changedCnt == 0 );
  }else {
    last;
  }
} # while(1);

unless ( exists $options{'quiet'} ){
  print "# Runtime: " . elapsedTime('main') . "\n";
}
  
&cleanup( $outdir, $conFile );
exit;

sub saveNewCons {
  my $conFile = shift;
  my $consRecs = shift;

  while ( -e "$conFile.$backupIdx" ) {
    $backupIdx++;
  }
  rename "$conFile", "$conFile.$backupIdx";
  open OUT, ">$conFile" or die "Could not open up $conFile for writing!\n";
  foreach my $consID ( keys(%{$consRecs}) ) {
    my $leftHPad = "H"x($consRecs->{$consID}->{'leftHPad'});
    my $rightHPad = "H"x($consRecs->{$consID}->{'rightHPad'});
    print OUT ">$consID";
    if ( $consRecs->{$consID}->{'class'} ne "" ) {
      print OUT "#" . $consRecs->{$consID}->{'class'};
    }
    if ( $consRecs->{$consID}->{'desc'} ne "" ) { 
      print OUT "  " . $consRecs->{$consID}->{'desc'};
    }
    print OUT "\n" . $leftHPad . $consRecs->{$consID}->{'seq'} . $rightHPad . "\n"; 
  }
  close OUT;
}

sub cleanup {
  my $outdir = shift;
  my $conFile = shift;
  # cleanup
  foreach my $ext ( "nog", "nsg", "nsi", "nhr", "nin", "nsq", "nsd" ){
    unlink "$outdir/$conFile.$ext" if ( -s "$outdir/$conFile.$ext" );
  }
}

sub scoretotal {
  my $file = shift;

  my $cmScoreTot = 0;
  my $alignedBases = 0;
  my %uniqIDs = ();
  open SCIN,"<$file" or die;
  while (<SCIN>) {
    if ( /^\s*(\d+)\s+\d+\.\d+\s+\d+\.\d+\s+\d\.\d+\s+(\S+)\s+(\d+)\s+(\d+)/ ) {
      $cmScoreTot += $1;
      $uniqIDs{$2}++;
      $alignedBases += $4-$3+1;
    }
  }
  close SCIN;
  print "Unique aligned sequences: " . scalar( keys(%uniqIDs) ) . "\n";
  print "Total Crossmatch Score: $cmScoreTot\n";
  print "Per Base Average: " . sprintf( "%0.2f", ( $cmScoreTot / $alignedBases )) . "\n";
  return (scalar(keys(%uniqIDs)),$cmScoreTot,sprintf( "%0.2f", ( $cmScoreTot / $alignedBases )));
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

  my $consSeq = "";
  my $refSeq = "";
  my $Hleft = 0;
  my $Hright = 0;
  my @consbit = ();
  my $consline = "";
  my $consbit = "";
  my $diffStr = "";
  while (<ALI>) {
    # consensus     1 TCTTCTGATT----GGTTGGTGGTGAGGTAA-------CA-G...
    if (/^consensus\s+\d+\s+(\S+)/) {
      $consSeq .= $1; 
      $consbit = $1;  # Save for comparison of sequences
      $consline = $_;
      chomp $consline;
    # ref:rep       1 HHHHHTGATT----GGTTGGTGGTGAGGTAA-------CA-G...
    } elsif (/^(ref\:\S+\s+\d+\s+)(\S+)/) {
      $refSeq .= $2; 
      my $spacelen = length $1;
      my $refbit = $2;
      if ( $refbit =~ /^(H+)[ACGTNRYMKSW-]/ ) {
        $Hleft = length($1);
      }
      # Could span multiple stanzas
      if ( $Hleft > 0 && $refbit =~ /[^H]+(H+)\s*$/ ) {
        $Hright += length($1);
      }
      if ($refbit ne $consbit) {
        if ( $is5ExtDone && /^ref\:\S+\s+(\d+)\s+(H+)[ACGTN]/ ) {
          # We were done with 5'extension previously.  Any alignment
          # to H's in this round are purely for anchoring purposes and
          # should not contribute to the diff output
          $Hleft = length $2;
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

  # Process new consensus into parts
  my $consHLeft;
  my $consHRight;
  my $consCore;
  my $cTmp = $consSeq;
  $cTmp =~ s/-//g;
  if ( $cTmp =~ /^(\w{$Hright})(\w+)(\w{$Hleft})$/ ) {
    $consHLeft = $1;
    $consCore = $2;
    $consHRight = $3;
  }else {
    print STDERR "ERROR: Could not parse cons $consSeq into parts!  Hright = $Hright Hleft = $Hleft\n";
    exit(1);
  }
 
  return( $consHLeft, $consCore, $consHRight, $consSeq, $refSeq, $diffStr );
}



#
# Read in the consensus file and build up a data structure for use
# by the main code.
#
sub processConFile {
  my $conFile = shift;
  my $optHpadCnt = shift;

  my($filename, $dirs, $suffix) = fileparse($conFile);
  open CON,"<$conFile" or die "\n\nCould not open consensus file \'$conFile\' for reading!\n\n";
  my $numBuffers = 0;
  my $idx = 1;
  my %consRecs = ();
  my $id;
  my $seq;
  my $desc;
  my $class;
  my $seqsWithHpads = 0;
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
        if ( $optHpadCnt ) {
          $leftHPad = $optHpadCnt;
          $rightHPad = $optHpadCnt;
        }
        $seqsWithHpads++ if ( $leftHPad > 0 || $rightHPad > 0 );
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
    if ( $optHpadCnt ) {
      $leftHPad = $optHpadCnt;
      $rightHPad = $optHpadCnt;
    }
    $seqsWithHpads++ if ( $leftHPad > 0 || $rightHPad > 0 );
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

  return( $numRefineableCons, $numBuffers, \%consRecs, $seqsWithHpads );
}

##-------------------------------------------------------------------------##
## Use: my $string = elapsedTime( $index );
##
## Great little utility for measuring the elapsed
## time between one or more sections of perl code.
##
##-------------------------------------------------------------------------##
sub elapsedTime {
  my ( $TimeHistIdx ) = @_;
  if ( defined $TimeBefore{$TimeHistIdx} ) {
    my $DiffTime = time - $TimeBefore{$TimeHistIdx};
    $TimeBefore{$TimeHistIdx} = time;
    my $Min = int( $DiffTime / 60 );
    $DiffTime -= $Min * 60;
    my $Hours = int( $Min / 60 );
    $Min -= $Hours * 60;
    my $Sec = $DiffTime;
    my $timeStr = sprintf( "%02d:%02d:%02d", $Hours, $Min, $Sec );
    return "$timeStr (hh:mm:ss)";
  }
  else {
    $TimeBefore{$TimeHistIdx} = time;
    return 0;
  }
}



1;
