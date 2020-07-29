#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) alignAndCallConsensus.pl
##  Author:
##      Arian Smit         asmit@systemsbiology.org
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A script to align a consensus to a set of putative
##      family elements and recall the consensus based on the
##      linup of all alignments.
##

=head1 NAME

alignAndCallConsensus.pl - Align TE consensus to a set of elements and call consensus

=head1 SYNOPSIS

  alignAndCallConsensus.pl [-c(onsensus) <*.fa>] [-e(lements) <*.fa>] 
                 [-d(ivergencemax) #] [-sc(ore) #] [-mi(nmatch) #]
                 [-h(tml)] [-st(ockholm)] [-q(uiet)]
                 [-b(andwidth) #] [-cr(ossmatch)] [-rm(blast)]
                 [-re(fine)]

  Example:

    ./alignAndCallConsensus.pl

        Assumes that a consensus file named "rep" and an elements file 
        named "repseq" ( representative sequences ) exists in the 
        current working directory. 

  or

    ./alignAndCallConsensus.pl -c mycons.fa -e myinstances.fa


  This tool was original written by Arian Smit and known affectionately
  as "dothemsimple.pl". As part of a large set of small command line tools
  this script performs an alignment of a putative TE family consensus to a 
  representative set of elements for the family.  The alignments are then
  used to generate a transitively induced multiple alignment ( aka "line up)
  for the purposes of consensus calling.  In a sense this is a hill climbing
  optimisation.  The output consensus may be better than the input consensus
  but it still may not be optimal or even stable.  By applying this method
  iteratively the consensus will most often stabilize at a local ( and perhaps
  global ) maximum.

  The tool was written in a way to allow for multiple use cases:

  Consensus Extension:
     Provided a consensus from a de-novo TE discovery tool, this tool may
     be used to manually extend a consensus sequence.  The "H" IUB code
     has been repurposed for this task.  By prefixing/suffxing "H" characters
     to the input consensus, flanking bases may be pulled into the final 
     alignment. The matrices used by this utility score any match to an
     "H" positively (+3).  This has the effect of extending the alignment
     by the number of "H" on either end of the alignment.  This assumes that
     the -elements file contains flanking sequence in advance.

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

  OUTPUT:

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

Default off; given n Zs on an end of the reference, the program wil extend 
the consensus by n bases and pad the new consensus with n Zs.
-f 5 : 5' extension done (the 5' end won't grow nor get Zs attached)
-f 3 : 3' extension done
-f b : both extensions done

=item -h(tml)       

Generate an HTML visualization of the alignment using viewMSA.pl. 
NOTE: This can be a slow process for large alignments and can be
easily be called outside this script.

=item -st(ockholm)     

Creates a Stockholm format (.stk) file on top of the normal Linup alignment file

=item -rm(blast)

Use RMBlast instead of cross_match

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
my $phrapDir = "/usr/local/bin/phrap";
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
                    '-finishedext|f=s',
                    '-matrix|ma=s',
                    '-minmatch|mi=s',
                    '-score|sc=s',
                    '-stockholm|st',
                    '-outdir=s',
                    '-html|h',
                    '-quiet|q',
                    '-refine|re',
                    '-bandwidth|b=s',
                    '-crossmatch|cr',
                    '-rmblast|rm',
                    '-onlyBestAlignment',
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


# Output directory
my $outdir = cwd;

my $conFile = "rep";
if ( exists $options{'consensus'} ) {
  if ( -s $options{'consensus'} ) {
    $conFile = $options{'consensus'}; 
    my($filename, $dirs, $suffix) = fileparse($options{'consensus'});
    $outdir = $dirs;
  }else {
    print "\n\nConsensus file \"$conFile\" is missing or empty!\n\n";
    exit(1);
  }
}elsif ( ! -s $conFile )
{
  print "\n\nDefault consensus file \"rep\" is missing or empty! See -c option for details.\n\n";
  exit(1);
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
  print "\n\nDefault elements file \"rep\" is missing or empty! See -c option for details.\n\n";
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

##
## Main
##
unless ( $options{'quiet'} ) {
  print "##\n";
  print "## alignAndCallConsensus (aka dothemsimple)\n";
  print "##\n";
  print "# Engine: $engine  Matrix: $matSpec, Bandwidth: $bandwidth,\n";
  print "#                  Minmatch: $minmatch, Minscore: $minscore,\n";
  print "#                  Maxdiv: $maxdiv, GapInit: $gapInit,\n";
  print "#                  InsGapExt: $insGapExt, DelGapExt: $delGapExt\n";
  if ( $options{'refine'} ) {
    print "#                  Refine: true\n";
  }else {
    print "#                  Refine: false\n";
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
  } else
  {
    # Sort the results collection
    $resultCollection->sort(
      sub ($$) {
        $_[ 1 ]->getScore() <=> $_[ 0 ]->getScore();
       }
    );
 
    open OUT,">$outdir/out" or die "could not open 'out' for writing\n";
    my $removedMinDiv = 10000;
    my $removedMaxDiv = 0;
    my $removedDivCount = 0;
    my $removedDupCount = 0;
    my %queryIDUsed = ();
    for ( my $k = 0 ; $k < $resultCollection->size() ; $k++ ) {
      my $resultRef = $resultCollection->get( $k );
      if ( $resultRef->getPctDiverge() <= $maxdiv ) {
        if ( $options{'onlyBestAlignment'} )
        {
          if ( exists $queryIDUsed{$resultRef->getQueryName()} ) {
            $removedDupCount++;
            next;
          }
          $queryIDUsed{$resultRef->getQueryName()}++;
          print OUT "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
        }else {
          print OUT "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
        }
      }else {
        $removedDivCount++;
        $removedMinDiv = $resultRef->getPctDiverge() if ( $resultRef->getPctDiverge() < $removedMinDiv);
        $removedMaxDiv = $resultRef->getPctDiverge() if ( $resultRef->getPctDiverge() > $removedMaxDiv);
      }
    }
    if ( $removedDivCount > 0 && ! $options{'quiet'} ){
      print "# Removed $removedDivCount high divegence alignments ranging from: $removedMinDiv - $removedMaxDiv divergence\n";
    }
    if ( $removedDupCount > 0 && ! $options{'quiet'} ){
      print "# Removed $removedDupCount alignments to that match the same repseq instance\n";
    }
    close OUT;
  }
  
  # Creating new ali file
  system "$FindBin::RealBin/Linup -i $outdir/out > $outdir/ali";
  
  if ( $options{'html'} ){
    system "$FindBin::RealBin/viewMSA.pl -order start -malign $outdir/out.malign";
    #system "cp MultipleAlignment.html ~/public_html";
  }
  
  if ( $options{'stockholm'} ) {
    system "$FindBin::RealBin/Linup -i $outdir/out ~/Matrices/linupmatrix -stockholm > $outdir/ali.stk";
  }
  
  unless ( $options{'quiet'} ) {
  print "-----\n";
  system "$FindBin::RealBin/checkfordifsinalis $outdir/ali";
  system "$FindBin::RealBin/scoretotal.pl $outdir/out";
  system "$FindBin::RealBin/CntSubst -cg $outdir/out | grep kimura";
  print "-----\n";
  }
  
  # 
  # Determine if H-padding was used and how long the pads were
  #
  my $original = `cat $conFile`;
  $original =~ s/^(>.+\n)//;
  my $faline = $1;
  $original =~ s/[\s]+//g; # should eliminate the line break at the end too                    
  my $h1 = "";
  my $h2 = "";
  if ($original =~ s/^([hH]+)(\S+?)([hH]*)$/$2/ ) {
    $h1 = $1;
    $h2 = $3 if $3;
  } elsif ($original =~ s/^(\S+?)([hH]*)$/$1/ ) {
    $h2 = $2 if $2;
  }

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
  
  #not all Hs may have been aligned to, so calculate length on each site                       
  my ($newcons,$Hleft,$Hright) = ();
  my (@consbit, $consline) = ();
  my $consbit;
  open (IN, "$outdir/ali") or die;
  while (<IN>) {
    if (/^consensus\s+\d+\s+(\S+)/) {
      $newcons .= $1;
      $consbit = $1;
      $newcons =~ tr/-//d;
      $consline = $_;
      chomp $consline;
    } elsif (/^(ref\:\S+\s+\d+\s+)(\S+)/) {
      my $spacelen = length $1;
      my $refbit = $2;
      if ($refbit ne $consbit) {
        if ( $ext5done && /^ref\:\S+\s+(\d+)\s+(H+)/ ) {
          if ($1 <= length $h1) { # so we're not at the end of the seq
            $Hleft = length $2;
            $newcons =~ s/^\w{$Hleft}//;
          }
          my ($tempcons,$tempref) = ($consbit,$refbit);
          $tempcons =~ s/^\w{$Hleft}//;
          $tempref =~ s/^Z+//;
          next if $tempcons eq $tempref;
        } elsif ($ext3done) {
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
        #print "$consline\n$midline\n$_\n";
      }
    }
  }
  
  $newcons =~ s/(\w){$Hright}$// if $Hright && $ext3done;
  if ($newcons eq $original) {
    print "Consensus unchanged\n" unless ( $options{'quiet'} );
    last;
  } else {
    my $backupIdx = 1;
    while ( -e "$conFile.$backupIdx" && $backupIdx < 25 ) {
      $backupIdx++;
    }
    if ( $backupIdx == 25 ) {
      die "The number of saved consensus files ( $conFile ) has exceded 25!  Please remove old backups before proceeding\n";
    }
    rename "$conFile", "$conFile.$backupIdx";
    open OUT, ">$conFile" or die;
    print OUT "$faline$h1\n$newcons\n$h2\n";
    close OUT;
  }
  last if ( ! $options{'refine'} );
  $iterations++;
  if ( $iterations >= 5 ) {
    print "Consensus still changing after 5 iterations...giving up.\n" unless ( $options{'quiet'} );
    last;
  }
}
  
# cleanup
foreach my $ext ( "nog", "nsg", "nsi", "nhr", "nin", "nsq", "nsd" ){
  unlink "$outdir/$conFile.$ext" if ( -s "$outdir/$conFile.$ext" );
}

1;
