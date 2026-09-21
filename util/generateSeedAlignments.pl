#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) generateSeedAlignments
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Generate seed alignments using RepeatMasker alignments
##      of a consensus library against an assembly.
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

generateSeedAlignments - Generate a seed alignments from RM *.align output

=head1 SYNOPSIS

 generateSeedAlignments [-families "<id1> <id2> .."] [-consensusRF]
                        [-outSTKFile <*.stk>] [-taxon <ncbi_taxonomy_name>]
                        [[-consensi <*.fa> [-outTable <*.tsv>] [-outAlign <*.align>]]
                        [-parasail] [-assemblyID <id>] [-minAlignedLength #]
                        [-targetCopies #] [-minDepth #] [-threads #]
                        [-verbose][-noColor] [-prefixAssembly][-filterDFDecoys]
                        -assemblyFile <*.2bit>
                        <RepeatMasker *.align File>

=head1 DESCRIPTION

Reverse engineer seed alignments using the following pipline:

  o [optional] : Screen initial consensus library for
    possible short period tandem repeat families.  These
    sometimes find their way in de-novo produced output.

  o Run RepeatMasker with a consensus library against the 
    assembly from which the consensus library was derived.
    (NOTE: use the -a option to produce an alignment file )

  o Run this script using the RepeatMasker alignment output
    and the assembly in 2bit format to produce seed alignments
    for each family ( or a particular set ) in Stockholm format.
    NOTE: The RF line in each Stockholm file represents the 
    match states as defined by the consensus used by RepeatMasker.
    This means that it may not match what would be derived by 
    a consensus call on the MSA.  As such we use the "X/." symbols
    in the RF line to indicate to Dfam that either "-use_ref_pos"
    needs to be set or the RF line needs to be updated.
    NOTE: The script samples instance depth at one position
    every 10 bp of the consensus.  If any sampled position has
    fewer than -minDepth instances, the script adds a "Coverage" line to
    the "#=GF **" lines of that seed alignment.  The line gives
    the number of sampled positions with no instances
    ( uncovered ), the number below the minimum depth ( low,
    which counts the uncovered ones too ), the number sampled,
    the minimum depth, and the lowest and highest depth.

  o [optional] : If there is a high level of fragmentation in
    the original consensus library, run an extension algorithm
    ( RepeatAfterMe/ExtendAlign, or Arian's dothemsimple method ).

  o Refine the consensus...using alignAndCallConsensus.pl.  The
    sample of instances chosen by this script probably don't match
    the sample used to generate the original consensus fed to 
    RepeatMasker.  It is essential to take the sample sequences and
    iterate a consensus building process to simultaneously refine
    the alignment and the consensus until they stabilize.

  o Compress the alignment.  We store seed alignments in Dfam using
    A3M.......document this!



The options are:

=over 4

=item -families "<id1> <id2> .."

Only analyze a specific set of families from the RepeatMasker alignment file.
Give the identifier exactly as the alignment file carries it, including any
"#class" suffix, as in "AluY#SINE/Alu".

=item DEPRECATED: -nucleotideRF replaced with -consensusRF

=item -consensusRF

The RF line produced by this tool is by-default derived from the sequence
used by RepeatMasker to identify copies. The use of this new flag changes
this behaviour by instead using the consensus derived directly from the 
MSA itself. By the Dfam convention we use the "x/." symbols whenever the RF
line is not a true consensus of the MSA and the consensus residues otherwise.
NOTE: this script calls the consensus once, from the alignment as
RepeatMasker produced it.  It does not iterate between calling a consensus
and realigning the instances to it, so the RF line reflects a single pass.
Run alignAndCallConsensus.pl on the seed alignment to iterate.

=item -outSTKFile <*.stk>

Concatenate all seed alignments generated into one file.  If not specified the
default is to create individual Stockholm files for each family.

=item -consensi <*.fa>

The consensus library that was given to RepeatMasker, in FASTA format.  The
identifiers may carry a "#class" suffix.

Each seed alignment carries a "ConsCmp" line recording how far the consensus
called from the alignment has moved from the consensus RepeatMasker used.
Supply this option to make that comparison over the full length of the family
("src=library").  Without it the script falls back to the reference sequence
of the multiple alignment ("src=alignedRef"), which is rebuilt from the
instances and so covers only the consensus positions they aligned to.

"sub" is the number of positions where both consensi have a base and the
bases differ.  "amb" is the number where the called base is not A, C, G or
T; those are left out of "sub".  With -consensi, a position where the
RepeatMasker base is an ambiguity code is not counted as "sub" either.
"ins" is the number of bases only the called consensus has, and "del" the
number only the RepeatMasker consensus has.  "id" is the percentage of
positions where both have a base that are neither "sub" nor "amb".  "lens"
holds two lengths, the called consensus then the RepeatMasker consensus.
"change" is the first length minus the second, so "change=+12" means the
called consensus is 12 bp longer.

With this option the script also adds a "ConsCAF" line: the alignment of the
two consensi as one record in RepeatMasker's CAF format ( see
SearchResult.pm, whose parseFromCAF() reads it back ).  The RepeatMasker
consensus is the query, named "OLD", and the called consensus the subject,
named "NEW".  In the last field
"G/C" is an OLD G opposite a NEW C, "+..+" encloses bases only NEW has, and
"-..-" encloses bases only OLD has.

This option also enables -outTable and -outAlign.

=item -outTable <*.tsv>

Write one tab separated row per family comparing the library consensus with
the consensus called from the seed alignment: family name, CpG counts before
and after, N counts before and after, lengths before and after, the number of
unambiguous substitutions, and the substitution, deletion and insertion
percentages.  Requires -consensi.

=item -outAlign <*.align>

Write the global alignment of each library consensus to the consensus called
from the seed alignment, in cross_match format.  Requires -consensi.

=item -parasail

Align each library consensus to the called consensus with parasail_aligner
( https://github.com/jeffdaily/parasail ) instead of the built-in Perl
aligner.  The Perl aligner's time and memory grow with the product of the two
consensus lengths.  In one test a 2.7 kb family took 34 seconds and 4.7 GB
with the Perl aligner and under a second with parasail_aligner.  The script
looks for parasail_aligner in $PARASAIL_DIR/bin and $PARASAIL_DIR if that
environment variable is set, and on the PATH otherwise.  It exits with an
error if it does not find the program.

Both aligners use the same matrix and gap penalties.  Where several
alignments tie for the best score they can choose differently, which moves a
gap along a run of repeated bases in the -outAlign output.  The ConsCmp and
-outTable counts can change too, though they did not in a test of five
families.  The aligners also differ when the alignment starts with a gap of
more than about 2 kb, as it can when one consensus is much longer than the
other.  The Perl aligner charges too little for that gap and
parasail_aligner does not.

=item -prefixAssembly

Prefix every instance identifier in the Stockholm output with the assembly
name, giving "<assembly>:<sequence>:<start>-<end>_<orient>".  Use it when
combining seed alignments from more than one assembly.

=item -taxon <ncbi_taxonomy_name>

Default taxon for to set in the Stockholm files for each family.

=item -assemblyID <id>

The identifier to attach to the family instance ranges.  If not specified
the default is to use the name of the assembly file.

=item -verbose

Include many more details in the log output.

=item -assemblyFile <*.2bit>

A two bit file containing the assembly that was RepeatMasked.  The
script checks the sequence of every alignment against it and drops
those that do not match.

=item -minAlignedLength

The minimum size a repeat instance must be to include in a seed alignment.
[Default = 30]

=item -targetCopies #

The number of copies to choose for each family.  The script ranks the copies
of a family by Kimura divergence: the least diverged three quarters first,
longest consensus span first, then the most diverged quarter in the same
order.  It takes copies from the first group until it has this many.  After
that, and for the most diverged quarter throughout, it takes a copy only if
the copy covers a sampled consensus position still below -minDepth.  A
family can therefore end up with more copies than this.
[Default = 500]

=item -minDepth #

The number of copies that should cover each sampled consensus position.  The
script samples one position every 10 bp of the consensus.  It keeps adding
copies that cover a position below this depth until every sampled position
reaches it or the copies run out.  For a family that falls short the script
prints a warning in the log and adds a "Coverage" line to its seed alignment.
[Default = 10]

=item -threads #

The number of families to build at once.  The script reads the alignment file
in one process whatever this is set to.  With more than one it then builds
each family in a child process.  A child starts out sharing the parent's copy
of the alignments, but much of it does not stay shared.  In a whole genome
run with a 13.6 GB parent the three children sampled reached 4 to 5 GB of
their own while building.  Reference count updates on the objects a child
reads are one likely cause.  The causes have not been measured.  With one
thread the
output files and the log list the families in sorted order.  With more than
one they list them in the order they finish.  Either way the script adds
each family to the output files as soon as it is built.  With more than one,
if a family fails the script lets the families already running finish, adds
them to the output files, and exits with an error.  More than one with -consensi
requires -parasail, because the Perl aligner's memory ( 4.7 GB for one
2.7 kb family ) is too much to run several at once.
[Default = 1]

=item -filterDFDecoys

In some cases it is desirable to include Dfam families in the search in order to
screen out matches to the existing family.  This option filters out alignemnts
to these families and doesn't generate seed alignments for anything with the
ID like "DF####..."

=item -noColor

Do not use escape codes to color the coverage depth log output.

=item -version

Print the version of this program and exit.

=back

=head1 ALSO

RepeatMasker, RepeatModeler, Dfam.org

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
#TODO:
#
#  - What primarily influences the lack of coverage in families?
#	* Is it that when competed through RepeatMasker some
#          related families steal instances?  If so, why doesn't 
#          RepeatModeler's competition with previous round consensi
#          reduce this?
#  - This script should refine the seed alignment/consensus before
#    creating a stockholm file.
#  - Store the divergence of the family as determined from these whole
#    genome runs.
#  - Make sure the stockholm contains the minimal fields for import into Dfam
#  - Update Classification Translation Table

#
# Module Dependence
#
use strict;
use sort 'stable';
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use RepModelConfig;
use MultAln;
use EMBL;
use SeedAlignment;
use SeedAlignmentCollection;
use SequenceSimilarityMatrix;             
use NeedlemanWunschGotohAlgorithm;
use File::Temp qw/ tempfile tempdir /;
use IO::Handle;
use POSIX ();
use File::Basename;
use Time::HiRes qw( gettimeofday tv_interval);

#
# RepeatMasker Dependencies
#
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use SearchResult;
use SearchResultCollection;
use CrossmatchSearchEngine;

#
# Paths
#
my $ucscToolsDir = $RepModelConfig::configuration->{'UCSCTOOLS_DIR'}->{'value'};

#
# Version
#
my $Version = $RepModelConfig::VERSION;

my %TimeBefore = ();

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',                # print out the version and exit
                    '-verbose',
                    '-families=s',
                    '-assemblyFile=s',
                    '-assemblyID=s',
                    '-consensi=s',
                    '-outTable=s',
                    '-outAlign=s',
                    '-taxon=s',
                    '-consensusRF',
                    '-outSTKFile=s',
                    '-prefixAssembly',
                    '-filterDFDecoys',
                    '-noColor',
                    '-parasail',
                    '-targetCopies=i',
                    '-minDepth=i',
                    '-threads=i',
                    '-minAlignedLength=i'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) )
{
  usage();
}

sub usage
{
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} )
{
  print "$Version\n";
  exit;
}

# 
# Experimental: Require seed alignment instances to have a particular set of diagnostic
#               sites. E.g for Aluyk3:
#  my @diagSites = ( [248, "C"], [252, "G"], [263, "A"], [236, "A"], [57, "C"], [86, "T"], [96, "C"] );
my @diagSites = ();

my %onlyTheseFamilies = ();
if ( defined $options{'families'} )
{
  foreach my $value ( split(/\s+/,$options{'families'} ) ) {
    $onlyTheseFamilies{$value} = 1;
  }
}

if ( ! $options{'assemblyFile'} ) {
  print "\nMust supply an assemblyFile parameter!\n\n";
  usage();
}

if ( $options{'outSTKFile'} ) {
  if ( -e $options{'outSTKFile'} ) {
    die "\nOutput file already exists $options{'outSTKFile'}.  Please remove and re-run program\n\n";
  }
}

my %consensi = ();
if ( $options{'consensi'} ) {
  open IN,"<$options{'consensi'}" or die "Could not open consensi file $options{'consensi'} for reading!\n";
  my $seq = "";
  my $id = "";
  while (<IN>) {
    if ( /^>(\S+)/ ) {
      my $tmpID = $1;
      # Family IDs are looked up without the "#class" suffix
      $tmpID =~ s/#.*//;
      if ( $seq ) {
        $consensi{$id} = $seq;
      }
      $id = $tmpID;
      $seq = "";
      next;
    }  
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  if ( $seq ) {
    $consensi{$id} = $seq;
  }
  close IN;
}

if ( ! $options{'consensi'} ) {
  print "\nWARNING: -consensi was not supplied.  The consensus comparison\n"
      . "         recorded for each family ( the ConsCmp line ) will fall back to the\n"
      . "         reference sequence of the multiple alignment, which is rebuilt from\n"
      . "         the instances and therefore covers only the consensus positions they\n"
      . "         aligned to.  Supply the library given to RepeatMasker with -consensi\n"
      . "         to compare over the full length of each consensus.\n\n";
}

my $alignFile = $ARGV[0];
if ( ! -s $alignFile ) {
  print "\nError: Missing RepeatMasker alignment file!\n\n";
  usage();
}

if (($options{'outTable'} || $options{'outAlign'}) && ! $options{'consensi'} ) {
  print "\nError: Options -outTable and -outAlign require the use of the -consensi option!\n\n";
  usage();
}

my $tableFH;
if ( $options{'outTable'} ) {
  open $tableFH,">$options{'outTable'}" or die "Could not open table file $options{'outTable'} for writing!\n";
}
my $alignFH;
if ( $options{'outAlign'} ) {
  open $alignFH,">$options{'outAlign'}" or die "Could not open align file $options{'outAlign'} for writing!\n";
}
my $NWMatrix = SequenceSimilarityMatrix->new();
$NWMatrix->parseFromFile( "$FindBin::Bin/../Matrices/linupmatrix" );

# With -parasail, find parasail_aligner under $PARASAIL_DIR if that is
# set, otherwise on the PATH.
my $parasailPrgm;
if ( $options{'parasail'} ) {
  my @candidates = ();
  if ( defined $ENV{'PARASAIL_DIR'} && $ENV{'PARASAIL_DIR'} ne "" ) {
    @candidates = ( "$ENV{'PARASAIL_DIR'}/bin/parasail_aligner",
                    "$ENV{'PARASAIL_DIR'}/parasail_aligner" );
  }else {
    @candidates = map { "$_/parasail_aligner" } grep { $_ ne "" } split( /:/, $ENV{'PATH'} );
  }
  ( $parasailPrgm ) = grep { -f $_ && -x $_ } @candidates;
  if ( ! $parasailPrgm ) {
    die "\nError: -parasail was given but parasail_aligner was not found "
        . ( $ENV{'PARASAIL_DIR'} ? "under PARASAIL_DIR ( $ENV{'PARASAIL_DIR'} )"
                                 : "on the PATH.  Set PARASAIL_DIR to the parasail install directory" )
        . "\n\n";
  }
  if ( $options{'consensi'} ) {
    print "Consensus comparison aligner: $parasailPrgm\n";
  }else {
    print "\nWARNING: -parasail has no effect without -consensi.\n\n";
  }
}

# The minimum sequence length to include in the seed alignment ( in bp ).
my $minAlignedLength = 30;
$minAlignedLength = $options{'minAlignedLength'} if ( $options{'minAlignedLength'} );

# The number of copies to choose for each family, and the depth every
# sampled consensus position should reach.
my $targetSampleCount = 500;
my $minDepth          = 10;
foreach my $opt ( 'targetCopies', 'minDepth' ) {
  if ( defined $options{$opt} && $options{$opt} < 1 ) {
    print "\nError: -$opt must be 1 or greater!\n\n";
    usage();
  }
}
$targetSampleCount = $options{'targetCopies'} if ( defined $options{'targetCopies'} );
$minDepth          = $options{'minDepth'}     if ( defined $options{'minDepth'} );

# The number of families to build at once
my $threads = 1;
if ( defined $options{'threads'} ) {
  if ( $options{'threads'} < 1 ) {
    print "\nError: -threads must be 1 or greater!\n\n";
    usage();
  }
  $threads = $options{'threads'};
}
# The Perl aligner took 4.7 GB on one 2.7 kb family, which is too much
# to run several of at once.
if ( $threads > 1 && $options{'consensi'} && ! $options{'parasail'} ) {
  print "\nError: -threads greater than 1 with -consensi requires -parasail!\n\n";
  usage();
}

elapsedTime("full_program");

#
# Parse the alignment file
#
elapsedTime("reading");
#
# The parser hands each alignment to this callback and stores nothing
# itself.  Alignments to Simple/Low/short families, Dfam decoys and
# families outside the -families set are dropped here before they take
# memory.  On whole-genome input that is most of the file.
#
my $resultCollection = SearchResultCollection->new();
my $numFilteredAtParse = 0;
my $keepResult = sub {
  my $result = shift;
  my $familyName = $result->getSubjName();
  if ( $familyName =~ /\#Simple|\#Low|short/ ) {
    $numFilteredAtParse++;
    return;
  }
  if ( $options{'filterDFDecoys'} && $familyName =~ /^DF\d\d\d\d\d.*/ ) {
    $numFilteredAtParse++;
    return;
  }
  if ( %onlyTheseFamilies && ! exists $onlyTheseFamilies{$familyName} ) {
    $numFilteredAtParse++;
    return;
  }
  $resultCollection->add( $result );
};
my $ALIGN;
if ( $alignFile =~ /\.gz$/ ) {
  open $ALIGN,"gunzip -c $alignFile|" or die "Could not run gunzip on $alignFile: $!\n";
}else {
  open $ALIGN,"<$alignFile" or die "Could not open $alignFile for reading: $!\n";
}
CrossmatchSearchEngine::parseOutput( searchOutput => $ALIGN, callback => $keepResult );
close $ALIGN;
print "Alignment file read in: " . elapsedTime("reading") . "\n";
print "  " . $resultCollection->size() . " alignments kept, $numFilteredAtParse skipped at parse time\n";


#
# Scrub the alignments
#
#   * Remove simple/tandem/low_complexity and "short" alignments
#
#   * Calculate the consensus size for each family and warn if 
#     inconsistent in the *.align file.  The most frequently reported
#     size wins.
#
#   * Handle some alignment artefacts generated by RepeatMasker's
#     processing of search-engine data.
#
elapsedTime("scrubbing");
my %consSizeByID = ();
my %consSizeCounts = ();
my $numBadRMAlignData = 0;
my $numShort = 0;
my $numLongXStretch = 0;
my %invalid = ();
for ( my $i = 0 ; $i < $resultCollection->size() ; $i++ ) {
  my $result = $resultCollection->get( $i );
  my $familyName = $result->getSubjName();

  # Tally reported consensus sizes
  $consSizeCounts{$familyName}{ $result->getSubjEnd() + $result->getSubjRemaining() }++;

  #
  # Handle some strange cases with RM alignment data:
  #
  my $querySeq = $result->getQueryString();
  my $subjSeq  = $result->getSubjString();

  # Alignments that begin/end with a double gap aligned characters
  # e.g:
  #
  #    1 ---ACAAT 4
  #    2 ---ACAAT 5
  #
  # This is an old bug and probably doesn't show up in modern
  # RM datasets.  This is not recoverable...just skip them.
  if (    ( $querySeq =~ /^-/ && $subjSeq =~ /^-/ )
       || ( $subjSeq =~ /-$/ && $querySeq =~ /-$/ ) )
  {
    print "Did not anticipate double gaps:\n"
        . $result->toStringFormatted( SearchResult::AlignWithQuerySeq )
        . "\n";
    $numBadRMAlignData++;
    $invalid{$i} = 1;
    next;
  }

  # Alignments that start in a gap.  This is typical of artifically
  # broken up alignments.  This is recoverable.  E.g:
  #
  #   Query 1 ------ACC 3
  #   Subj  2 ACACCTACC 11
  #
  # Just chew back to be:
  #   
  #   Query 1 ACC 3
  #   Subj  8 ACC 11
  # 
  my ( $gapChars ) = ( $querySeq =~ /^(\-+).*/ );
  # Query Gap Start
  if ( $gapChars )
  {
    $querySeq = substr( $querySeq, length( $gapChars ) );
    $subjSeq  = substr( $subjSeq,  length( $gapChars ) );
    $result->setQueryString( $querySeq );
    $result->setSubjString( $subjSeq );
    if ( $result->getOrientation() eq "C" )
    {
      $result->setSubjEnd( $result->getSubjEnd() - length( $gapChars ) );
    }else {
      $result->setSubjStart( $result->getSubjStart() + length( $gapChars ) );
    }
  }
  # Subj Gap Start
  ( $gapChars ) = ( $subjSeq =~ /^(\-+).*/ );
  if ( $gapChars )
  {
    my $gapQSeq = substr( $querySeq, 0, length( $gapChars ) );
    $gapQSeq =~ s/X//g;
    $querySeq = substr( $querySeq, length( $gapChars ) );
    $subjSeq  = substr( $subjSeq,  length( $gapChars ) );
    $result->setQueryString( $querySeq );
    $result->setSubjString( $subjSeq );
    # NOTE: querySequence may contain X's ....watch out that we don't
    #       count these as query positions.
    $result->setQueryStart( $result->getQueryStart() + length( $gapQSeq ) );
  }
  ( $gapChars ) = ( $subjSeq =~ /(\-+)$/ );
  # Subj Gap End
  if ( $gapChars )
  {
    my $gapQSeq =
       substr( $querySeq, length( $querySeq ) - length( $gapChars ) );
    $gapQSeq =~ s/X//g;
    $querySeq =
        substr( $querySeq, 0, length( $querySeq ) - length( $gapChars ) );
    $subjSeq =
        substr( $subjSeq, 0, length( $subjSeq ) - length( $gapChars ) );
    $result->setQueryString( $querySeq );
    $result->setSubjString( $subjSeq );
    $result->setQueryEnd( $result->getQueryEnd() - length( $gapQSeq ) );

  }
  # Query Gap End
  ( $gapChars ) = ( $querySeq =~ /(\-+)$/ );
  if ( $gapChars )
  {
    $querySeq =
        substr( $querySeq, 0, length( $querySeq ) - length( $gapChars ) );
    $subjSeq =
        substr( $subjSeq, 0, length( $subjSeq ) - length( $gapChars ) );
    $result->setQueryString( $querySeq );
    $result->setSubjString( $subjSeq );
    if ( $result->getOrientation() eq "C" )
    {
      $result->setSubjStart( $result->getSubjStart() + length( $gapChars ) );
    }else {
      $result->setSubjEnd( $result->getSubjEnd() - length( $gapChars ) );
    }
  }

  # Alignments that have an internal repeat ( long string of Xs )
  # These are problematic.  If we include the internal repeat it will
  # get built into the family model as a mosaic.  For now just
  # report the number of times it occurs.
  if ( $querySeq =~ /X{10}/ ) {
    $numLongXStretch++; 
    $invalid{$i} = 1;
    next;
  }
  #if ( $querySeq =~ /X{10}/ && $options{'assembly'} )
  #{
  #  my $xStart     = $result->getQueryStart() - 1;
  #  my $twoBitFile = $options{'assembly'};
  #  my $chr         = $result->getQueryName();
  #  my $newQuerySeq = "";
  #  while ( $querySeq =~ /([^X]+)(X{10,}+)/ig )
  #  {
  #    my $prefixSeq    = $1;
  #    my $xSeq         = $2;
  #    my $prefixSeqLen = $prefixSeq;
  #    $prefixSeqLen =~ s/-//g;
  #    $xStart += length( $prefixSeqLen );
  #    my $xEnd = $xStart + length( $xSeq );
  #    $newQuerySeq .= $prefixSeq;
  #    my $replSeq = `$ucscToolsDir/twoBitToFa $twoBitFile:$chr:$xStart-$xEnd stdout`;
  #    $xStart += length( $xSeq );
  #    $replSeq =~ s/^>[^\n\r]+[\n\r]+//;
  #    $replSeq =~ s/[\n\r]//g;
  #    $newQuerySeq .= $replSeq;
  #  }
  #  $newQuerySeq .= substr( $querySeq, length( $newQuerySeq ) );
  #  $result->setQueryString( $newQuerySeq );
  #  $querySeq = $newQuerySeq;
  #}
  # Alignments that still contain an 'X' character
  if ( $querySeq =~ /X/ )
  {
    # This was due to a bug in 4.0.5 and earlier.  Shouldn't happen
    # after that.  Perhaps warn?
    $numBadRMAlignData++;
    $invalid{$i} = 1;
    next;
  }
    
  # Test for minimum length
  my $qs = $querySeq;
  $qs =~ s/-//g;
  if ( length( $qs ) < $minAlignedLength )
  {
    $numShort++;
    $invalid{$i} = 1;
    next;
  }
}

# Settle on one consensus size per family
foreach my $familyName ( keys %consSizeCounts ) {
  my $counts = $consSizeCounts{$familyName};
  my @sizes = sort { $counts->{$b} <=> $counts->{$a} || $a <=> $b } keys %$counts;
  $consSizeByID{$familyName} = $sizes[0];
  if ( @sizes > 1 ) {
    print "WARN: $familyName has more than one reported size: "
        . join( ", ", map { "$_ (x$counts->{$_})" } @sizes )
        . "; using $sizes[0]\n";
  }
}
undef %consSizeCounts;

print "$numBadRMAlignData bad RepeatMasker alignment data\n";
print "$numShort sequences were too short ( < $minAlignedLength bp )\n";
print "$numLongXStretch sequences with internal masked region (>= 10 Xs in a row).\n";
print "Scrubbing RM data in: " . elapsedTime("scrubbing") . "\n";

#
# Validate that the sequence coordinates match the given assembly file
#
elapsedTime("validating");
my @survivorIdx = grep { ! $invalid{$_} } 0 .. $resultCollection->size() - 1;
my @survivors   = map { $resultCollection->get( $_ ) } @survivorIdx;
my %failedRef   = map { ( $_ => 1 ) } validateAgainstAssembly( \@survivors );
my $numFailedSeqValidation = 0;
for ( my $j = 0 ; $j < @survivors ; $j++ ) {
  next unless ( $failedRef{ $survivors[$j] } );
  $invalid{ $survivorIdx[$j] } = 1;
  $numFailedSeqValidation++;
}
undef @survivors;
undef %failedRef;
print "$numFailedSeqValidation sequences failed validation against the assembly.\n";
print "Validating sequences in: " . elapsedTime("validating") . "\n";



# Generate alignment pointers organized by family name
my %alignByID = ();
for ( my $i = 0 ; $i < $resultCollection->size() ; $i++ ) {
  unless ( $invalid{$i} ) {
    my $result = $resultCollection->get( $i );
    my $familyName = $result->getSubjName();
    push @{$alignByID{$familyName}}, $result;
  }
}

undef $resultCollection;
print "  Total Families:  " . scalar( keys( %alignByID ) ) . " ( excluding simple/low )\n";

my $consLen              = 0;
my $totalAttemptedBuilds = 0;

# Set the identifier for the assembly to be used in Stockholm file.  If not provided
# explicitly use the filename of the assembly (sans path).
my $assemblyName = $options{'assemblyFile'};
$assemblyName = basename($assemblyName);
$assemblyName =~ s/\.2bit//i;
$assemblyName = $options{'assemblyID'} if ( $options{'assemblyID'} );

my $totalBuilt = 0;
my $noAlign = 0;
my $numNoCov = 0;
my $numPoorCov = 0;

##
## Build the seed alignment for one family.  With -threads greater than 1
## this runs in a child process, with stdout and the output files pointed
## at that family's own temporary files.
##
sub buildFamily
{
  my $id = shift;

  elapsedTime("family_build_time");
  $totalAttemptedBuilds++;
  my $countInGenome = scalar(@{$alignByID{$id}});
  $consLen = $consSizeByID{$id};
  print "Working on $id ( length=$consLen, $countInGenome in assembly )\n";
  elapsedTime("outlier_detection");

  # Sort alignments for this family by divergence (ascending) and find the median and
  # quartile values.  Sort keys are computed once per alignment and the sort is
  # stable, so ties keep their input order.
  my @sortedByDiv = map  { $_->[1] }
                    sort { $a->[0] <=> $b->[0] }
                    map  { [ $_->getPctKimuraDiverge(), $_ ] } @{$alignByID{$id}};
  my $n = scalar(@sortedByDiv);
  my $medianDiv;
  if ( $n % 2 )
  {
    # Odd
    $medianDiv = $sortedByDiv[int($n/2)]->getPctKimuraDiverge();
  }else {
    # Even
    $medianDiv = ($sortedByDiv[int($n/2)-1]->getPctKimuraDiverge() +
               $sortedByDiv[int($n/2)]->getPctKimuraDiverge()
              ) / 2;
  }
  my $quartileIdx = int( 3 * $n / 4 );
  my $quartileDiv = $sortedByDiv[$quartileIdx]->getPctKimuraDiverge();
  print "  Median Divergence = $medianDiv\n";
  print "  3rd Quartile Divergence = $quartileDiv\n";

  # Remove elements that are in top divergence quartile
  my @outliers = splice(@sortedByDiv, $quartileIdx);

  # Sort both lists by consensus span (descending)
  my $bySpan = sub {
    return map  { $_->[1] }
           sort { $b->[0] <=> $a->[0] }
           map  { [ $_->getSubjEnd() - $_->getSubjStart(), $_ ] } @_;
  };
  @outliers = $bySpan->( @outliers );
  print "  " . scalar(@outliers) . " outliers were moved to the end of the priority list\n";
  my @data = $bySpan->( @sortedByDiv );
  my $outlierStartIdx = scalar(@data);

  # Precedence:
  #     Elements within first 3 quartiles of kimura divergence
  #     Long elements
  #     Elements that cover a low-sample-depth region
  push @data, @outliers;

  # Coverage is sampled at one position every 10 bp of the consensus
  # ( positions 1, 11, 21, ... ).  A copy covers a sampled position when
  # the position lies within its consensus range.
  my $numBins      = int( ( $consLen + 9 ) / 10 );
  my @sampledDepth = ( 0 ) x $numBins;
  my $lowBins      = $numBins;   # sampled positions still below $minDepth
  my $idx          = 0;
  my $sampleCount  = 0;
  my @chosen       = ();
  my %seen         = ();
  my $numBadReportedConsLen = 0;
  my $numDiagMismatch = 0;
  my $dupsFound  = 0;
  print "   - Sorting and data preparation : " . elapsedTime("outlier_detection") . "\n";
  elapsedTime("instance_selection");

  # Walk the priority list until the sample budget is met and every
  # sampled position has reached $minDepth.
  while ( @data )
  {
    if ( $sampleCount >= $targetSampleCount && $lowBins == 0 )
    {
      print "  Coverage reached!\n";
      last;
    }
    my $result = shift @data;
    $idx++;

    # This was warned about in the scrubbing routine...is this necessary?
    my $calcCLen = $result->getSubjEnd() + $result->getSubjRemaining();
    if ( $calcCLen != $consLen )
    {
      $numBadReportedConsLen++;
    }

    # Experimental ( currently not implemented for multi-family use )
    # This only keeps elements that contain the correct diagnostic
    # sites and bases.
    if ( @diagSites )
    {
      my $qry = $result->getQueryString();
      my $sbj = $result->getSubjString();

      my @sPosToQBase = ();
      my $sIdx = $result->getSubjStart(); # one based
      $sIdx = $result->getSubjEnd() if ( $result->getOrientation() eq "C" );
      for ( my $i = 0; $i < length($sbj); $i++ )
      {
        my $qBase = substr($qry, $i, 1);
        if ( $qBase eq "-" )
        {
          next;
        }
        if ( $result->getOrientation() eq "C" )
        {
          $qBase =~ tr/ACGT/TGCA/;
          $sPosToQBase[$sIdx] = $qBase;
          $sIdx--;
        }else {
          $sPosToQBase[$sIdx] = $qBase;
          $sIdx++;
        }
      }

      my $failed = 0;
      foreach my $site ( @diagSites )
      {
        if ( $sPosToQBase[$site->[0]] ne $site->[1] )
        {
          $failed = 1;
        }
      }

      if ( $failed )
      {
        $numDiagMismatch++;
        next;
      }
    }

    # Skip exact duplicates of a range already seen.  Fragments of one
    # copy have different ranges and are kept.
    my $key =   $result->getQueryName()  . ":"
              . $result->getQueryStart() . "-"
              . $result->getQueryEnd();
    if ( $seen{$key} ) {
      $dupsFound++;
      next;
    }
    $seen{$key}++;

    # Sampled positions covered by this copy: the first at or after its
    # consensus start and the last at or before its consensus end.
    my $iStart = int( ( $result->getSubjStart() + 8 ) / 10 );
    my $iEnd   = int( ( $result->getSubjEnd() - 1 ) / 10 );
    $iEnd = $numBins - 1 if ( $iEnd > $numBins - 1 );
    my $addIt = 0;
    for ( my $i = $iStart ; $i <= $iEnd ; $i++ )
    {
      if ( $sampledDepth[$i] < $minDepth )
      {
        $addIt = 1;
        last;
      }
    }
    if ( ( $idx <= $outlierStartIdx && $sampleCount < $targetSampleCount ) || $addIt )
    {
      for ( my $i = $iStart ; $i <= $iEnd ; $i++ )
      {
        $sampledDepth[$i]++;
        $lowBins-- if ( $sampledDepth[$i] == $minDepth );
      }
      # RepeatMasker writes masked bases as X.  Use the IUPAC N instead.
      my $str = $result->getQueryString();
      $str =~ s/X/N/g;
      $result->setQueryString( $str );
      $str = $result->getSubjString();
      $str =~ s/X/N/g;
      $result->setSubjString( $str );
      push @chosen, $result;
      $sampleCount++;
    }
  }
  print "   - Selecting instances : " . elapsedTime( "instance_selection" ) . "\n";

  my $noCovExamples = 0;
  my $minCovDepth = 10000000000;
  my $maxCovDepth = 0;
  my $idxStrLen = length($consLen);
  print "  ";
  for ( my $i = 0 ; $i < $numBins ; $i++ )
  {
    my $idxPos = ( $i * 10 ) + 1;
    if ( $options{'noColor'} ) {
      print "[" . sprintf("%$idxStrLen"."s",$idxPos) . "]=" . sprintf("%5s", $sampledDepth[$i]) . ", ";
    }else {
      if ( $sampledDepth[$i] >= $minDepth ) {
        print "[" . sprintf("%$idxStrLen"."s",$idxPos) . "]=" . sprintf("%5s", $sampledDepth[$i]) . ", ";
      }elsif (  $sampledDepth[$i] > 0 ) {
        # Yellow
        print "[" . sprintf("%$idxStrLen"."s",$idxPos) . "]=" . 
              "\033[33m" . sprintf("%5s", $sampledDepth[$i]) . "\033[0m" . ", ";
      }else {
        # Red
        print "[" . sprintf("%$idxStrLen"."s",$idxPos) . "]=" . 
              "\033[31m" . sprintf("%5s", $sampledDepth[$i]) . "\033[0m" . ", ";
      }
    }
    if ( ($i+1) % 10 == 0 ) {
      print "\n  ";
    }
    $minCovDepth = $sampledDepth[$i] if ( $sampledDepth[$i] < $minCovDepth ); 
    $maxCovDepth = $sampledDepth[$i] if ( $sampledDepth[$i] > $maxCovDepth ); 
    $noCovExamples++ if ( $sampledDepth[$i] < 1 );
  }
  print "\n";

  print "  Stats:\n";
  print "    Coverage depth range: $minCovDepth to $maxCovDepth ( from sampled positions )\n";
  print "    $sampleCount of $countInGenome were chosen for the multiple alignment\n";
  print "    $numBadReportedConsLen had differing data about consensus length ( RM artifact )\n";
  print "    $dupsFound duplicate alignments (same exact range)\n";
  if ( @diagSites ) {
    print "    $numDiagMismatch had mismatches to the specified diagnostic sites\n";
  }

  if ( $noCovExamples )
  {
    print "  *** Some regions are not covered! ***\n";
    $numNoCov++;
  }

  if ( $minCovDepth < $minDepth )
  {
    print "  *** Some regions did not reach the min coverage depth of $minDepth  ***\n";
    $numPoorCov++;
  }

  if ( $sampleCount == 0 ) {
    print "WARNING: $id is being skipped because there are no alignments?!??\n";
    $noAlign++;
    return;
  }

  my $resultCol = SearchResultCollection->new();
  foreach my $result ( @chosen )
  {
    if ( $options{'prefixAssembly'} ) {
      $result->setQueryName( $assemblyName . ":" . $result->getQueryName() );
    }
    $resultCol->add( $result );
  }

  my $mAlign = MultAln->new( searchCollection          => $resultCol,
                             searchCollectionReference => MultAln::Subject );

  #
  # Compare the consensus called from this alignment with the consensus
  # RepeatMasker used to find the instances.  The consensus is called once
  # here, with no iteration between calling a consensus and realigning the
  # instances to it.  Where the two consensi disagree, the alignment no
  # longer supports the consensus it was built from.
  #
  # With -consensi the comparison is against the library sequence over its
  # full length.  Without it the only copy of the RepeatMasker consensus on
  # hand is the reference sequence of the multiple alignment, which is
  # reconstructed from the instances and therefore covers just the consensus
  # positions they aligned to.  The two are reported under different labels
  # because they do not measure the same thing.
  #
  my $calledCons = $mAlign->consensus();
  my $libID = $id;
  $libID = $1 if ( $id =~ /(\S+)#.*/ );
  my $libCons = $consensi{$libID};
  my $ungappedCons = $calledCons;
  $ungappedCons =~ s/[- ]//g;
  my %consCmp = ();
  my @libCmp  = ();
  if ( defined $libCons && $libCons ne "" ) {
    my %opts = ( familyName => $libID, oldSeq => uc( $libCons ),
                 newSeq => uc( $ungappedCons ), SSMatrixObj => $NWMatrix );
    $opts{'alignFH'} = $alignFH if ( $alignFH );
    $opts{'parasailPrgm'} = $parasailPrgm if ( $parasailPrgm );
    @libCmp = compareConsensi( %opts );
    $consCmp{'source'}      = "library";
    $consCmp{'substituted'} = $libCmp[ 6 ];
    $consCmp{'inserted'}    = $libCmp[ 10 ];
    $consCmp{'deleted'}     = $libCmp[ 11 ];
    $consCmp{'ambiguous'}   = $libCmp[ 12 ];
    $consCmp{'aligned'}     = $libCmp[ 13 ];
    $consCmp{'caf'}         = $libCmp[ 14 ];
    $consCmp{'refLen'}      = $libCmp[ 4 ];
  }else {
    print "WARNING: $libID is not in the consensus library given with -consensi.  "
        . "Comparing against the aligned consensus positions instead.\n"
        if ( $options{'consensi'} );
    %consCmp = compareConsToReference( $mAlign->getReferenceSeq(), $calledCons );
    $consCmp{'source'} = "alignedRef";
    $consCmp{'refLen'} = $consLen;
  }
  my $consPctId = "NA";
  $consPctId = sprintf( "%0.2f%%",
        100 * ( $consCmp{'aligned'} - $consCmp{'substituted'} - $consCmp{'ambiguous'} )
            / $consCmp{'aligned'} )
      if ( $consCmp{'aligned'} );
  print "  Called consensus vs RepeatMasker consensus ( $consCmp{'source'} ):\n";
  print "    $consCmp{'substituted'} substitutions, $consCmp{'inserted'} inserted, "
      . "$consCmp{'deleted'} deleted, $consCmp{'ambiguous'} ambiguous "
      . "over $consCmp{'aligned'} compared positions ( $consPctId identity )\n";
  print "    $consCmp{'aligned'} of $consCmp{'refLen'} consensus positions were compared\n"
      if ( $consCmp{'refLen'} );
  # Few instances give low identity on their own, since the consensus is
  # then just those instances.  Low identity flags a build for review only
  # where the family has many instances.

  if ( ( $options{'outTable'} || $options{'outAlign'} ) && @libCmp ) {
    my ($oldCpGCount, $newCpGCount, $oldNCount, $newNCount, $len_old, $len_new,
        $nonAmbigSubCount, $pctSub, $pctDel, $pctIns) = @libCmp;
    print "    CpG = $oldCpGCount -> $newCpGCount, N = $oldNCount -> $newNCount,\n";
    print "    Length = $len_old -> $len_new (" . ($len_new - $len_old) . "), NonAmbig Substititions = $nonAmbigSubCount,\n";                
    print "    Sub/Del/Ins = $pctSub $pctDel $pctIns\n";                 
    if ( $tableFH ) {
      print $tableFH "$libID\t$oldCpGCount\t$newCpGCount\t$oldNCount\t$newNCount\t$len_old\t$len_new\t$nonAmbigSubCount\t$pctSub\t$pctDel\t$pctIns\n";
    }
  }

  my $sanitizedID = $id;
  $sanitizedID =~ s/[\(\)]//g; 
  my $class = "Unknown";
  if ( $sanitizedID =~ /(.*)\#(.*)/ )
  {
    $sanitizedID = $1;
    $class = $2;
  }
  # handle names that have "/" in them
  $sanitizedID =~ s/\//_/g;
  
  $totalBuilt++;
  if ( $options{'consensusRF'} ) {
    $mAlign->toSTK(
                      filename         => "$sanitizedID.stk",
                      consRF => 1,
                      id               => $sanitizedID
      );
  } else
  {
    $mAlign->toSTK(
                    filename => "$sanitizedID.stk",
                    id       => $sanitizedID
    );
  }

  # Patch up the stockholm
  my $stockholmFile = SeedAlignmentCollection->new();
  open my $IN, "<$sanitizedID.stk"
      or die "Could not open up stockholm file $sanitizedID.stk for reading!\n";
  $stockholmFile->read_stockholm( $IN );
  close $IN;
  unlink("$sanitizedID.stk");

  my $seedAlign   = $stockholmFile->get( 0 );
  $seedAlign->setClassification(&RMClassToDfam($class));
  if ( $options{'taxon'} ) {
    $seedAlign->addClade($options{'taxon'});
  }
  my $desc = "Seed alignments generated from RepeatMasker annotations using generateSeedAlignments.pl. ".
             "The median Kimura divergence for the family is $medianDiv, $sampleCount were chosen from $countInGenome identified in " . 
             "the $assemblyName assembly.";
  $seedAlign->setComments($desc);
  # "lens" is the called consensus length then the RepeatMasker consensus
  # length, and "change" the first minus the second.
  my $calledLen = length( $ungappedCons );
  my $lenDiff   = $calledLen - $consCmp{'refLen'};
  $lenDiff = "+$lenDiff" if ( $lenDiff > 0 );
  $desc = "Source:gsa, mDiv=$medianDiv, $assemblyName:$countInGenome\n"
        . "ConsCmp: src=$consCmp{'source'}, sub=$consCmp{'substituted'}, "
        . "ins=$consCmp{'inserted'}, del=$consCmp{'deleted'}, "
        . "amb=$consCmp{'ambiguous'}, "
        . "lens=$calledLen/$consCmp{'refLen'}, change=$lenDiff, "
        . "id=$consPctId";
  # Only the comparison against the library has a pairwise alignment
  $desc .= "\nConsCAF: $consCmp{'caf'}" if ( defined $consCmp{'caf'} );
  # Repeat the log's coverage warnings in the Stockholm file, which is
  # read without the log.
  if ( $noCovExamples || $minCovDepth < $minDepth )
  {
    $desc .= "\nCoverage: uncovered=$noCovExamples, low=$lowBins, "
           . "sampled=$numBins, minDepth=$minDepth, "
           . "depth=$minCovDepth-$maxCovDepth";
  }
  $seedAlign->setCuratorComments($desc);

  if ( $options{'outSTKFile'} ) 
  {
    open OUT,">>".$options{'outSTKFile'} or die;
  }else { 
    open OUT,">$sanitizedID.stk" or die;
  }

  print OUT "" . $seedAlign->toString();
  close OUT;
  print "   - total build time : " . elapsedTime("family_build_time") . "\n";
}

# The families to build.  With one thread the output files list them in
# this order.
my @familyIDs = sort grep { ! $options{'families'} || exists $onlyTheseFamilies{$_} }
                     keys( %alignByID );

if ( $threads <= 1 ) {
  foreach my $id ( @familyIDs ) {
    buildFamily( $id );
  }
}else {
  #
  # Build up to $threads families at once, each in a forked child.  A
  # child starts with the parent's memory shared, not copied, so the
  # alignments read above are not duplicated up front.  Much of it does
  # not stay shared: in a whole genome run with a 13.6 GB parent the
  # children sampled reached 4 to 5 GB of their own while building.  A
  # child writes its
  # log, its part of each output file and its counters to files in
  # $pieceDir.  As each child finishes the parent prints its log and
  # appends its parts to the output files, so the families come out in
  # the order they finish.
  #
  my $pieceDir = tempdir( "tmpGenSeedsXXXXXXXX", DIR => ".", CLEANUP => 1 );
  my %pieceIdx = ();
  @pieceIdx{ @familyIDs } = ( 0 .. $#familyIDs );

  # Start the families with the most copies x consensus length first, so
  # that one long family is not left running alone at the end.
  my %cost  = map { $_ => scalar( @{ $alignByID{$_} } ) * $consSizeByID{$_} } @familyIDs;
  my @queue = sort { $cost{$b} <=> $cost{$a} || $a cmp $b } @familyIDs;

  my %running = ();
  my @failed  = ();
  while ( ( @queue && ! @failed ) || %running ) {
    while ( @queue && ! @failed && scalar( keys( %running ) ) < $threads ) {
      my $id    = shift @queue;
      my $piece = "$pieceDir/$pieceIdx{$id}";

      # Anything still buffered would be written again by the child
      STDOUT->flush();
      $tableFH->flush() if ( $tableFH );
      $alignFH->flush() if ( $alignFH );

      my $pid = fork();
      die "\nERROR: could not fork: $!\n" if ( ! defined $pid );
      if ( $pid == 0 ) {
        my $built = eval {
          open( STDOUT, ">", "$piece.log" ) or die "Could not open $piece.log: $!\n";
          $options{'outSTKFile'} = "$piece.stk" if ( $options{'outSTKFile'} );
          if ( $tableFH ) {
            open( $tableFH, ">", "$piece.tsv" ) or die "Could not open $piece.tsv: $!\n";
          }
          if ( $alignFH ) {
            open( $alignFH, ">", "$piece.align" ) or die "Could not open $piece.align: $!\n";
          }
          ( $totalAttemptedBuilds, $totalBuilt, $noAlign, $numNoCov, $numPoorCov ) = ( 0 ) x 5;
          buildFamily( $id );
          close $tableFH if ( $tableFH );
          close $alignFH if ( $alignFH );
          open my $CNT, ">", "$piece.cnt" or die "Could not open $piece.cnt: $!\n";
          print $CNT "$totalAttemptedBuilds $totalBuilt $noAlign $numNoCov $numPoorCov\n";
          close $CNT;
          close STDOUT;
          1;
        };
        print STDERR "\nERROR: building $id failed: $@\n" if ( ! $built );
        # Leave without Perl's normal teardown.  A plain exit frees every
        # object the child inherited, which writes to those pages and
        # gives the child its own copy of all the alignments.  In a test
        # each child of a 528 MB parent reached about 490 MB that way.
        # With _exit the children building small families stayed at 32
        # to 66 MB.
        STDOUT->flush();
        STDERR->flush();
        POSIX::_exit( $built ? 0 : 1 );
      }
      $running{$pid} = $id;
    }

    my $pid = wait();
    last if ( $pid == -1 );
    next if ( ! exists $running{$pid} );
    my $status = $?;
    my $id     = delete $running{$pid};
    my $piece  = "$pieceDir/$pieceIdx{$id}";
    appendPiece( "$piece.log", \*STDOUT );
    if ( $status != 0 || ! -s "$piece.cnt" ) {
      push @failed, $id;
      next;
    }

    open my $CNT, "<", "$piece.cnt" or die "Could not open $piece.cnt: $!\n";
    my @counts = split( " ", <$CNT> );
    close $CNT;
    unlink( "$piece.cnt" );
    $totalAttemptedBuilds += $counts[ 0 ];
    $totalBuilt           += $counts[ 1 ];
    $noAlign              += $counts[ 2 ];
    $numNoCov             += $counts[ 3 ];
    $numPoorCov           += $counts[ 4 ];

    if ( $options{'outSTKFile'} && -e "$piece.stk" ) {
      open my $STK, ">>", $options{'outSTKFile'}
          or die "Could not open $options{'outSTKFile'} for appending: $!\n";
      appendPiece( "$piece.stk", $STK );
      close $STK;
    }
    appendPiece( "$piece.tsv",   $tableFH ) if ( $tableFH );
    appendPiece( "$piece.align", $alignFH ) if ( $alignFH );
  }
  if ( @failed ) {
    die "\nERROR: the build failed for: " . join( ", ", @failed ) . "\n"
        . "  The output files are incomplete.\n\n";
  }
}

##
## Copy a child's output file onto the end of an open filehandle and
## remove the file.  A skipped family has no .stk file.
##
sub appendPiece {
  my ( $file, $toFH ) = @_;
  return if ( ! -e $file );
  open my $PIECE, "<", $file or die "Could not open $file for reading: $!\n";
  while ( <$PIECE> ) {
    print $toFH $_;
  }
  close $PIECE;
  unlink( $file );
}

print "\n\n";
print "Total Seeds alignments built: $totalBuilt out of $totalAttemptedBuilds\n";
print "    - Number with poor coverage areas ( depth < $minDepth ): $numPoorCov\n";
print "    - Number with no coverage areas: $numNoCov\n";
print "    - Number without any alignments: $noAlign\n";
print "Total runtime : " . elapsedTime("full_program") . "\n";
print "\n\n";

# All done
exit;

############################################################################################

##-------------------------------------------------------------------------##
## Use: my %counts = compareConsToReference( $referenceSeq, $calledCons );
##
##   Compare the consensus called from a multiple alignment with the
##   reference sequence of that alignment.  The reference is rebuilt from
##   the instances, so it holds only the consensus positions they aligned
##   to.  Use this when the library given to RepeatMasker is not available;
##   with -consensi, compare against the library sequence instead.  Both
##   sequences are in the same column coordinates, so the comparison is
##   column by column.  Returns a hash with these counts:
##
##     substituted  both called a base and the bases differ
##     matched      both called the same base
##     inserted     the reference has a gap and the alignment called a base
##     deleted      the reference has a base and the alignment called a gap
##     ambiguous    the alignment called something other than A, C, G or T
##     aligned      columns where both called a base ( matched + substituted
##                  + ambiguous )
##     uncovered    columns with no reference sequence to compare against
##-------------------------------------------------------------------------##
sub compareConsToReference {
  my ( $refSeq, $consSeq ) = @_;

  my %counts = ( substituted => 0, matched   => 0, inserted  => 0,
                 deleted     => 0, ambiguous => 0, aligned   => 0,
                 uncovered   => 0 );

  # The called consensus may run past the end of the reference when an
  # instance extends beyond the consensus RepeatMasker used.
  my $len = length( $refSeq );
  $len = length( $consSeq ) if ( length( $consSeq ) > $len );
  $refSeq  .= " " x ( $len - length( $refSeq ) );
  $consSeq .= " " x ( $len - length( $consSeq ) );

  for ( my $i = 0 ; $i < $len ; $i++ ) {
    my $r = uc( substr( $refSeq,  $i, 1 ) );
    my $c = uc( substr( $consSeq, $i, 1 ) );
    my $rIsBase = ( $r ne "-" && $r ne " " );
    my $cIsBase = ( $c ne "-" && $c ne " " );
    if ( ! $rIsBase && $r ne "-" ) {
      # No reference sequence to compare against
      $counts{'uncovered'}++;
    }elsif ( $rIsBase && $cIsBase ) {
      $counts{'aligned'}++;
      if ( $c !~ /[ACGT]/ ) {
        $counts{'ambiguous'}++;
      }elsif ( $c ne $r ) {
        $counts{'substituted'}++;
      }else {
        $counts{'matched'}++;
      }
    }elsif ( $rIsBase ) {
      $counts{'deleted'}++;
    }elsif ( $cIsBase ) {
      $counts{'inserted'}++;
    }
  }
  return %counts;
}


##-------------------------------------------------------------------------##
## Use: my @failed = validateAgainstAssembly( \@results );
##
##   Extract the genomic range of each result from the assembly and
##   compare it to the aligned query sequence.  Returns the results
##   whose sequence does not match.  Dies if twoBitToFa fails, since
##   a partial run would let unchecked ranges through.
##-------------------------------------------------------------------------##
sub validateAgainstAssembly {
  my ( $results ) = @_;
  return () unless ( @$results );

  my ( $tfh, $tfilename ) = tempfile( "tmpGenSeedsXXXXXXXX", DIR => ".", UNLINK => 1 );
  my %expected = ();
  my %rangeOf  = ();
  foreach my $result ( @$results ) {
    my $range = $result->getQueryName() . ":"
              . ( $result->getQueryStart() - 1 ) . "-"
              . $result->getQueryEnd();
    # twoBitToFa requires at least four BED columns
    print $tfh $result->getQueryName() . "\t"
             . ( $result->getQueryStart() - 1 ) . "\t"
             . $result->getQueryEnd() . "\t$range\n";
    my $qs = $result->getQueryString();
    $qs =~ s/-//g;
    $expected{$range} = uc( $qs );
    $rangeOf{$result} = $range;
  }
  close $tfh;

  my %failed = ();
  my $cmd = "$ucscToolsDir/twoBitToFa -bedPos -bed=$tfilename $options{'assemblyFile'} stdout";
  open my $IN, "$cmd|" or die "Could not run $cmd: $!\n";
  my $id  = "";
  my $seq = "";
  my $check = sub {
    return if ( $id eq "" );
    if ( ! exists $expected{$id} ) {
      print "ERROR: twoBitToFa returned $id, which was not requested\n";
    }elsif ( $expected{$id} ne uc( $seq ) ) {
      $failed{$id} = length( $seq );
    }
  };
  while ( <$IN> ) {
    if ( /^>(\S+)/ ) {
      $check->();
      $id  = $1;
      $seq = "";
      next;
    }
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  $check->();
  close $IN
      or die "\nERROR: twoBitToFa exited with status " . ( $? >> 8 )
           . " while validating ranges against $options{'assemblyFile'}\n\n";
  unlink( $tfilename );

  my @failedResults = ();
  foreach my $result ( @$results ) {
    my $range = $rangeOf{$result};
    next unless ( exists $failed{$range} );
    push @failedResults, $result;
    if ( $options{'verbose'} ) {
      print "Invalid sequence $range ( aligned seq length = "
          . length( $expected{$range} )
          . ", assembly seq length = $failed{$range} )\n";
    }
  }
  return @failedResults;
}


sub compareConsensi {
  my %nameValueParams = @_;

  my $familyName = $nameValueParams{'familyName'};
  my $ss_matrix = $nameValueParams{'SSMatrixObj'};
  my $matrix = $nameValueParams{'MatrixObj'};
  my $oldSeq = $nameValueParams{'oldSeq'};
  my $newSeq = $nameValueParams{'newSeq'};

  my $searchResult;
  if ( $nameValueParams{'parasailPrgm'} && $oldSeq ne "" && $newSeq ne "" ) {
    $searchResult = parasailGlobalAlign(
                  program        => $nameValueParams{'parasailPrgm'},
                  querySeq       => $oldSeq,
                  subjectSeq     => $newSeq,
                  matrix         => $ss_matrix,
                  gapOpenPenalty => -25,
                  gapExtPenalty  => -5
                );
  }
  if ( ! defined $searchResult ) {
    $searchResult = NeedlemanWunschGotohAlgorithm::search(
                  querySeq   => $oldSeq,
                  subjectSeq => $newSeq,
                  matrix         => $ss_matrix,
                  insOpenPenalty => -25,
                  insExtPenalty  => -5,
                  delOpenPenalty => -25,
                  delExtPenalty  => -5
                );
  }

  # A count of zero comes back from s/// as the empty string
  my $oldCpGCount = () = ( $oldSeq =~ /CG/ig );
  my $newCpGCount = () = ( $newSeq =~ /CG/ig );
  my $oldNCount   = () = ( $oldSeq =~ /N/ig );
  my $newNCount   = () = ( $newSeq =~ /N/ig );

  my $qs = $searchResult->getQueryString();
  my $ss = $searchResult->getSubjString();
  my $sub     = 0;
  my $ins     = 0;
  my $del     = 0;
  my $ambig   = 0;
  my $aligned = 0;
  for ( my $i = 0; $i < length($qs); $i++ ) {
     my $qbase = uc(substr($qs,$i,1));
     my $sbase = uc(substr($ss,$i,1));
     # The query is the consensus RepeatMasker used, the subject the one
     # called from the seed alignment.
     if ( $qbase eq "-" ) {
       $ins++ if ( $sbase ne "-" );
       next;
     }
     if ( $sbase eq "-" ) {
       $del++;
       next;
     }
     $aligned++;
     if ( $sbase !~ /[ACGT]/ ) {
       $ambig++;
     }elsif ( $qbase =~ /[ACGT]/ && $sbase ne $qbase ) {
       $sub++;
     }
  }

  # The alignment as one CAF record ( see SearchResult::_toCAF ), with the
  # RepeatMasker consensus as the query, "OLD".
  $searchResult->setQueryName("OLD");
  $searchResult->setSubjName("NEW");
  my $caf = $searchResult->toStringFormatted( SearchResult::CompressedAlignFormat );
  $caf =~ s/[\n\r]+$//;

    if ( $nameValueParams{'alignFH'} ) {
    if ( ref( $nameValueParams{'alignFH'} ) =~ /GLOB|FileHandle|IO::File/ ) {
      my $FH = $nameValueParams{'alignFH'};
      print $FH "Global Alignment: $familyName\n";
      print $FH "" . $searchResult->toStringFormatted( SearchResult::AlignWithQuerySeq );
      print $FH "\n\n";
    }
  }

  return $oldCpGCount, $newCpGCount, $oldNCount, $newNCount, length($oldSeq), length($newSeq), $sub, $searchResult->getPctDiverge(), 
         $searchResult->getPctDelete(), $searchResult->getPctInsert(), $ins, $del, $ambig, $aligned, $caf;
}

##
## Globally align two sequences with parasail_aligner and return a
## SearchResult with the fields NeedlemanWunschGotohAlgorithm::search()
## fills in.
##
## The matrix and gap penalties are the same as the Perl aligner's.
## Where several alignments tie for the best score the two aligners
## may pick different ones, so gaps can sit at different positions.
##
## Returns undef if parasail_aligner reports no alignment.
##
sub parasailGlobalAlign {
  my %parameters = @_;

  my $prgm     = $parameters{'program'};
  my $querySeq = $parameters{'querySeq'};
  my $subjSeq  = $parameters{'subjectSeq'};
  my $matrix   = $parameters{'matrix'};
  # parasail takes gap penalties as positive numbers
  my $gapOpen  = -$parameters{'gapOpenPenalty'};
  my $gapExt   = -$parameters{'gapExtPenalty'};

  # File::Temp removes the directory when $tmpDirObj goes out of scope.
  # The tempdir() form waits for the program to exit normally, which a
  # child started by -threads never does.
  my $tmpDirObj = File::Temp->newdir( "tmpGenSeedsXXXXXXXX", DIR => ".", CLEANUP => 1 );
  my $tmpDir    = $tmpDirObj->dirname();

  # parasail scores a pair as matrix[ database base ][ query base ], and
  # the Perl aligner as matrixHash{ subject base, query base }.  So with
  # the subject as the database the rows go out unchanged.  parasail uses
  # the "*" column for a base the matrix lacks.  It is 0, which is what the
  # Perl aligner adds for such a base.
  my @alphabet = @{ $matrix->{'alphabetArray'} };
  open my $MAT, ">$tmpDir/matrix" or die "Could not open $tmpDir/matrix for writing!\n";
  print $MAT "  " . join( " ", map { sprintf( "%4s", $_ ) } @alphabet, "*" ) . "\n";
  foreach my $rowBase ( @alphabet ) {
    print $MAT "$rowBase " . join( " ", map { sprintf( "%4d", $_ ) }
                   ( map { $matrix->{'matrixHash'}->{ $rowBase, $_ } } @alphabet ), 0 ) . "\n";
  }
  print $MAT "* " . join( " ", map { sprintf( "%4d", 0 ) } @alphabet, "*" ) . "\n";
  close $MAT;

  open my $QRY, ">$tmpDir/query.fa" or die "Could not open $tmpDir/query.fa for writing!\n";
  print $QRY ">query\n$querySeq\n";
  close $QRY;
  open my $SBJ, ">$tmpDir/subject.fa" or die "Could not open $tmpDir/subject.fa for writing!\n";
  print $SBJ ">subject\n$subjSeq\n";
  close $SBJ;

  # -x turns off the filter that aligns only pairs sharing an exact match
  # ( 7 bases by default ).  Use the 32 bit kernel: the score of a
  # consensus several kb long does not fit in 16 bits.  parasail_aligner
  # counts stdin as a third input.  It refused to run with stdin inherited
  # from this script or redirected from /dev/null, and ran with it closed.
  my $cmd = "$prgm -x -t 1 -a nw_trace_striped_32 -m $tmpDir/matrix "
          . "-o $gapOpen -e $gapExt -q $tmpDir/query.fa -f $tmpDir/subject.fa "
          . "-O SAM -g $tmpDir/out.sam <&- 2>&1";
  my $cmdOutput = `$cmd`;
  die "\nERROR: parasail_aligner exited with status " . ( $? >> 8 ) . "\n"
      . "  command: $cmd\n  output: $cmdOutput\n" if ( $? );

  my ( $cigar, $score );
  open my $SAM, "<$tmpDir/out.sam" or die "Could not open $tmpDir/out.sam for reading!\n";
  while ( <$SAM> ) {
    next if ( /^\@/ );
    my @fields = split( /\t/ );
    # parasail_aligner 2.6.2 writes an alignment that scores exactly 0 as
    # an unmapped read, with no CIGAR.  The caller falls back to the Perl
    # aligner.
    if ( $fields[ 1 ] & 4 ) {
      close $SAM;
      return undef;
    }
    $cigar = $fields[ 5 ];
    $score = $1 if ( /\tAS:i:(-?\d+)/ );
    last;
  }
  close $SAM;
  die "\nERROR: could not parse the parasail_aligner output for command: $cmd\n"
      unless ( defined $score && defined $cigar && $cigar =~ /^(\d+[=XMID])+$/ );

  # In the SAM record the read is the query.  "I" is a query base with
  # no subject base, and "D" the reverse.
  my $queryAlignment = "";
  my $subjAlignment  = "";
  my $qPos = 0;
  my $sPos = 0;
  while ( $cigar =~ /(\d+)([=XMID])/g ) {
    my ( $len, $op ) = ( $1, $2 );
    if ( $op eq "I" ) {
      $queryAlignment .= substr( $querySeq, $qPos, $len );
      $subjAlignment  .= "-" x $len;
      $qPos += $len;
    }elsif ( $op eq "D" ) {
      $queryAlignment .= "-" x $len;
      $subjAlignment  .= substr( $subjSeq, $sPos, $len );
      $sPos += $len;
    }else {
      $queryAlignment .= substr( $querySeq, $qPos, $len );
      $subjAlignment  .= substr( $subjSeq, $sPos, $len );
      $qPos += $len;
      $sPos += $len;
    }
  }
  die "\nERROR: the parasail_aligner alignment covers $qPos of " . length( $querySeq )
      . " query bases and $sPos of " . length( $subjSeq ) . " subject bases: $cigar\n"
      if ( $qPos != length( $querySeq ) || $sPos != length( $subjSeq ) );

  # The same statistics NeedlemanWunschGotohAlgorithm::search() reports
  my $insCnt  = ( $queryAlignment =~ tr/-// );
  my $percIns = sprintf( "%0.1f", ( $insCnt * 100 ) / length( $subjSeq ) );
  my $delCnt  = ( $subjAlignment =~ tr/-// );
  my $percDel = sprintf( "%0.1f", ( $delCnt * 100 ) / length( $querySeq ) );
  my $sub = 0;
  for ( my $i = 0; $i < length( $queryAlignment ); $i++ ) {
    my $qbase = uc( substr( $queryAlignment, $i, 1 ) );
    my $sbase = uc( substr( $subjAlignment, $i, 1 ) );
    next if ( $qbase eq "-" || $sbase eq "-" );
    $sub++ if ( $qbase =~ /[ACGT]/ && $sbase ne $qbase );
  }
  my $percSub = sprintf( "%0.1f", ( $sub * 100 ) / length( $querySeq ) );

  return SearchResult->new(
                            queryName      => "query",
                            subjName       => "subject",
                            pctInsert      => $percIns,
                            pctDelete      => $percDel,
                            subjStart      => 1,
                            subjEnd        => length( $subjSeq ),
                            queryStart     => 1,
                            queryEnd       => length( $querySeq ),
                            score          => $score,
                            pctDiverge     => $percSub,
                            subjRemaining  => 0,
                            queryRemaining => 0,
                            orientation    => "",
                            queryString    => $queryAlignment,
                            subjString     => $subjAlignment
  );
}



sub RMClassToDfam {
  my $rmClass = lc(shift);

#
# NOTE: This is a temporary translation table from the RepeatMasker classication
#       scheme to the Dfam_consensus one.  It's temporary for two reasons.  First
#       the RM scheme has a one-one mapping with the Dfam_consensus scheme at this
#       stage but that is not guaranteed to last.  Second, we intend to use the
#       new scheme in the classifier at some point making it unnecessary to do this
#       mapping or maintain this table.
#
# RepeatMasker type/subtype to Dfam Classification Lineage
# Autogenerated from Dfam/Server/exportClassification.py on 2022-03-28 12:41:14.896216
my %rmToDfamClass = (
"other" => "Other",
"artefact" => "Artifact",
"segmental" => "Segmental_Duplication",
"low_complexity" => "Low_Complexity",
"other/dna_virus" => "Accidental;Normally_Non-integrating_Virus",
"unknown" => "Interspersed_Repeat;Unknown",
"simple_repeat" => "Tandem_Repeat;Simple",
"satellite" => "Tandem_Repeat;Satellite",
"unknown/centromeric" => "Interspersed_Repeat;Unknown;Centromeric",
"rna" => "Interspersed_Repeat;Pseudogene;RNA",
"satellite/centromeric" => "Tandem_Repeat;Satellite;Centromeric",
"satellite/macro" => "Tandem_Repeat;Satellite;Macro",
"satellite/y-chromosome" => "Tandem_Repeat;Satellite;Y-chromosomal",
"satellite/acromeric" => "Tandem_Repeat;Satellite;Acromeric",
"satellite/w-chromosome" => "Tandem_Repeat;Satellite;W-chromosomal",
"satellite/subtelomeric" => "Tandem_Repeat;Satellite;Subtelomeric",
"dna" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase",
"rc" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Helicase",
"line" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE",
"scrna" => "Interspersed_Repeat;Pseudogene;RNA;scRNA",
"rrna" => "Interspersed_Repeat;Pseudogene;RNA;rRNA",
"trna" => "Interspersed_Repeat;Pseudogene;RNA;tRNA",
"snrna" => "Interspersed_Repeat;Pseudogene;RNA;snRNA",
"dna/crypton" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton",
"dna/p" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;P_Element",
"dna/kolobok" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok",
"dna/hat" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT",
"dna/ginger" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Ginger",
"dna/zisupton" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Zisupton",
"dna/zator" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Zator",
"dna/pif" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger",
"dna/merlin" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Merlin",
"dna/mule" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like",
"dna/dada" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Dada",
"dna/novosib" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Novosib",
"dna/tcmar" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner",
"dna/is3eu" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;IS3EU",
"dna/maverick" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;DNA_Polymerase;Maverick",
"rc/helitron-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Helicase;Helitron-2",
"rc/helitron" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Helicase;Helitron-1",
"retroposon" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter",
"sine" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE",
"line/penelope" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements",
"ltr" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element",
"unknown/tate" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;Viper-group;TATE",
"dna/crypton-s" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-S",
"dna/crypton-r" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-R",
"dna/crypton-v" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-V",
"dna/crypton-f" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-F",
"dna/crypton-c" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-C",
"dna/crypton-i" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-I",
"dna/crypton-x" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-X",
"dna/crypton-a" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-A",
"dna/crypton-h" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-H",
"dna/p-fungi" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;P_Element;Fungi-specific_Branch",
"dna/kolobok-hydra" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;Hydra-specific_Branch",
"dna/kolobok-h" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;Kolobok-H",
"dna/kolobok-e" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;Kolobok-E",
"dna/kolobok-t2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;T2",
"dna/hat-restless" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Restless",
"dna/hat-blackjack" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Blackjack",
"dna/hat-hobo" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hobo",
"dna/hat-tag1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Tag1",
"dna/hat-hatw" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hATw",
"dna/hat-tip100" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Tip100",
"dna/hat-hat6" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT6",
"dna/hat-hat19" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT19",
"dna/hat-hat1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT1",
"dna/hat-hatx" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hATx",
"dna/hat-pegasus" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Pegasus",
"dna/hat-hatm" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hATm",
"dna/hat-hat5" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT5",
"dna/hat-ac" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Activator",
"dna/hat-charlie" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Charlie",
"dna/cmc-transib" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;Transib",
"dna/cmc" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC",
"dna/sola-1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Sola;Sola-1",
"dna/sola-3" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Sola;Sola-3",
"dna/sola-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Sola;Sola-2",
"dna/pif-isl2eu" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;ISL2EU",
"dna/pif-harbinger" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;Harbinger",
"dna/pif-spy" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;Spy",
"dna/pif-harbs" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;HarbS",
"dna/academ-h" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Academ;Academ-H",
"dna/academ-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Academ;Academ-2",
"dna/academ-1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Academ;Academ-1",
"dna/piggybac-a" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PiggyBac_Group;PiggyBac-A",
"dna/piggybac" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PiggyBac_Group;PiggyBac",
"dna/piggybac-x" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PiggyBac_Group;PiggyBac-X",
"dna/mule-f" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;F",
"dna/mule-mudr" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;MuDR",
"dna/mule-nof" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;NOF",
"dna/casposons" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;DNA_Polymerase;Casposon",
"dna/tcmar-tc2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Tc2",
"dna/tcmar-tc4" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc4",
"dna/tcmar-mogwai" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Mogwai",
"dna/tcmar-stowaway" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Stowaway",
"dna/tcmar-mariner" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Mariner",
"dna/tcmar-cweed" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Cweed",
"dna/tcmar-sagan" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Sagan",
"dna/tcmar-gizmo" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Gizmo",
"dna/tcmar-ant1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Ant1",
"dna/tcmar-tc1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc1",
"dna/tcmar-isrm11" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;ISRm11",
"dna/tcmar-m44" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;m44",
"retroposon/r4-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;R4-derived",
"retroposon/l1-dep" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-dependent",
"retroposon/l1-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-derived",
"retroposon/l2-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L2-derived",
"retroposon/rte-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;RTE-derived",
"retroposon/i-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;I-derived",
"sine/7sl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter",
"sine/5s" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter",
"ltr/trim" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;TRIM",
"ltr/pao" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Bel-Pao",
"ltr/copia" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Ty1-Copia",
"line/cre-odin" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;Odin",
"line/cre" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;CRE",
"line/cre-ambal" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;Ambal",
"line/genie" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Genie",
"dna/cmc-enspm" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;EnSpm",
"dna/cmc-mirage" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;Mirage",
"dna/tcmar-tigger" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Tigger",
"dna/tcmar-fot1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Fot1",
"dna/tcmar-pogo" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Pogo",
"retroposon/sva" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-dependent;SVA",
"ltr/cassandra" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae",
"ltr/gypsy" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Gypsy",
"ltr/dirs" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;DIRS",
"ltr/viper" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;Viper-group;Viper",
"ltr/ngaro" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;Ngaro",
"line/cre-2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;CRE;CRE-2",
"line/cre-1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;CRE;CRE-1",
"line/proto1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;Proto-1",
"line/proto2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;Proto-2",
"line/rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like",
"dna/cmc-chapaev" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;Chapaev_group;Chapaev",
"dna/cmc-chapaev-3" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;Chapaev_group;Chapaev-3",
"sine/trna-5s" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_and_5S_RNA;No_or_Unknown_Core;Unknown_LINE-dependent",
"sine/ceph" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;Unknown_Promoter;Ceph-core;RTE-end",
"sine/core-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;Unknown_Promoter;MIR-core;RTE-end",
"sine/core" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;Unknown_Promoter;MIR-core;Unknown_LINE-dependent",
"sine/trna-v-cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;CR1-end",
"sine/trna-v" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;Unknown_LINE-dependent",
"sine/trna-sauria-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;L2-end",
"sine/trna-sauria-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;RTE-end",
"sine/trna-sauria" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;Unknown_LINE-dependent",
"sine/trna-ceph-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Ceph-core;RTE-end",
"sine/trna-ceph" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Ceph-core;Unknown_LINE-dependent",
"sine/trna-l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;L1-dependent",
"sine/trna-i" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;I-end",
"sine/trna-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;L2-end",
"sine/trna-rex" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Rex-end",
"sine/trna-cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;CR1-end",
"sine/trna-jockey" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Jockey-end",
"sine/r1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;R1-end",
"sine/trna-tad1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Tad1_End",
"sine/trna-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;RTE-end",
"sine/trna-r2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;R2-end",
"sine/trna" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Unknown_LINE-dependent",
"sine/rte-bovb" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;BovB-end",
"sine/trna-deu-i" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;I-end",
"sine/trna-deu-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;L2-end",
"sine/trna-deu-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;RTE-end",
"sine/trna-deu-cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;CR1-end",
"sine/trna-deu" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;Unknown_LINE-dependent",
"sine/id" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No-core;L1-dependent",
"sine/trna-v-core-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V_and_MIR-core;L2-end",
"sine/trna-meta" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Meta-core;Unknown_LINE-dependent",
"sine/mir" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;L2-end",
"sine/trna-mermaid" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;Mermaid",
"sine/trna-core-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;RTE-end",
"sine/trna-core" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;Unknown_LINE-dependent",
"sine/5s-sauria-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Sauria-core;RTE-end",
"sine/5s-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;No_or_Unknown_Core;RTE-end",
"sine/5s-deu-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Deu-core;L2-end",
"sine/5s-deu" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Deu-core;Unknown_LINE-dependent",
"sine/5s-core-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;MIR-core;RTE-end",
"sine/u" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;U-RNA_Promoter;No_or_Unknown_Core;Unknown_LINE-dependent",
"sine/b4" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_and_7SL_RNA;No-core;L1-dependent",
"sine/trna-7sl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_and_7SL_RNA;No_or_Unknown_Core;Unknown_LINE-dependent",
"ltr/erv-foamy" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Spumaretrovirinae",
"ltr/erv" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae",
"ltr/caulimovirus" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Pararetroviridae;Caulimoviridae",
"line/r2-hero" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R2-like;Hero",
"line/r2-nesl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R2-like;NeSL",
"line/r2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R2-like;R2",
"line/l1-dre" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;DRE",
"line/l1-zorro" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;Zorro",
"line/dualen" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R4-like;Dualen",
"line/dong-r4" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R4-like;Dong-R4",
"line/deceiver" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R4-like;Deceiver",
"line/rte-orte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;ORTE",
"line/tad1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;Tad1",
"sine/alu" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter;No-core;L1-dependent;Alu",
"sine/b2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No-core;L1-dependent;B2",
"ltr/erv1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1",
"line/l1-tx1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;L1-group;Tx1",
"line/l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;L1-group;L1",
"line/rte-x" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;RTE-group;RTE-X",
"line/rte-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;RTE-group;RTE",
"line/rte-bovb" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;RTE-group;BovB",
"line/rex-babar" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;Rex-Babar",
"line/cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;CR1",
"ltr/erv-lenti" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;Lenti",
"ltr/erv4" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV4",
"ltr/ervk" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV2",
"ltr/ervl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV3",
"line/l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;L2-group;L2",
"line/cr1-zenon" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;CR1;Zenon",
"line/r1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;R1-subgroup;R1",
"line/r1-loa" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;R1-subgroup;LOA",
"line/i" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;I-group;I",
"line/i-jockey" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;I-group;Jockey",
"ltr/ervl-malr" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV3;MaLR",
"dna/mule-ricksha" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;Ricksha",
"ltr/dirs-q" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;DIRS;Q",
"dna/maverick-mavirus" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;DNA_Polymerase;Maverick;Maverick-Mavirus",
"dna/tcmar-is885" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;IS885",
"sine/trna-v-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;L2-end",
"sine/trna-v-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;RTE-end",
"dna/ginger-1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Ginger;Ginger-1",
"dna/ginger-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Ginger;Ginger-2",
"retroposon/sno" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;snoRNA",
"sine/erv1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1;SINE-like",
"sine/u-l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;U-RNA_Promoter;No_or_Unknown_Core;L1-dependent",
"retroposon/sno-l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;snoRNA;L1-derived",
"ple/chlamys" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Chlamys",
"ple/hydra" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Hydra",
"ple/naiad" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Naiad",
"ple/athena" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Athena",
"ple/poseidon" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Poseidon",
"ple/nematis" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Nematis",
"ple/neptune" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Neptune",
"ple/coprina" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Coprina",
);

  if ( exists $rmToDfamClass{$rmClass} ){
    return $rmToDfamClass{$rmClass};
  }else {
    return $rmToDfamClass{"unknown"};
  }
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
    my $DiffTime = tv_interval( $TimeBefore{$TimeHistIdx} );
    $TimeBefore{$TimeHistIdx} = [ gettimeofday() ];
    my $Min = int( $DiffTime / 60 );
    $DiffTime -= $Min * 60;
    my $Hours = int( $Min / 60 );
    $Min -= $Hours * 60;
    my $timeStr = sprintf( "%02d:%02d:%06.3f", $Hours, $Min, $DiffTime );
    return "$timeStr (hh:mm:ss.sss) Elapsed Time";
  }
  else {
    $TimeBefore{$TimeHistIdx} = [ gettimeofday() ];
    return 0;
  }
}

1;
