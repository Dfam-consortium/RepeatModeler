#!/usr/local/bin/perl 
##---------------------------------------------------------------------------##
##  File:
##      @(#) generateSeedAlignments
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
#
## ChangeLog
#
#    $Log$
#
###############################################################################
#
#

=head1 NAME

generateSeedAlignments - Generate a seed alignment from RM *.align output

TODO:
   - Option to generate a single vs multiple stockholm files
   - Store the divergence of the family as determined from these whole
     genome runs.
   - Make sure the stockholm contains the minimal fields for import into Dfam

=head1 SYNOPSIS

generateSeedAlignments [-element <id>] [-includeRef]
                       [-outSTKFile <*.stk>] [-taxon <ncbi_taxonomy_name]
                       [-assemblyID <id>] -assembly <*.2bit>
                       <RepeatMasker *.align File>

=head1 DESCRIPTION

The options are:

=over 4

=item -element <id>

Only analyze one family from the provided RepeatMasker alignment file.

=item -includeRef

blah blah blah

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
use File::Temp qw/ tempfile tempdir /;
use Time::HiRes qw( gettimeofday tv_interval);

#
# RepeatMasker Dependencies
#
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use SearchResult;
use SearchResultCollection;
use CrossmatchSearchEngine;

#
#
#
my $ucscToolsDir = "/home/rh105648e/ucscTools";

#
# Version
#
my $Version = "2.0";

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
                    '-element=s',
                    '-assembly=s',
                    '-taxon=s',
                    '-update',
                    '-includeRef',
                    '-outSTKFile=s',
                    '-includeRefWithGaps',
                    '-excludeExistingDfam=s'
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

# Playing around with developing a seed alignment for AluYk3 -- success!
#my @diagSites = ( [248, "C"], [252, "G"], [263, "A"], [236, "A"], [57, "C"], [86, "T"], [96, "C"] );
my @diagSites = ();

my $theseElements = "";
if ( defined $options{'element'} )
{
  $theseElements = $options{'element'};
}

if ( $options{'outSTKFile'} ) {
  if ( -e $options{'outSTKFile'} ) {
    die "Output file already exists $options{'outSTKFile'}.  Please remove and re-run program\n";
  }
}

my $alignFile = $ARGV[0];
if ( ! -s $alignFile ) {
  die "Error: Missing RepeatMasker alignment file!\n";
}

elapsedTime("reading");
my $resultCollection;
my $ALIGN;
if ( $alignFile =~ /.*\.gz$/ ) {
  open $ALIGN,"gunzip -c $alignFile|" or die;
  $resultCollection =
      CrossmatchSearchEngine::parseOutput( searchOutput => $ALIGN );
  close $ALIGN;
}else {
  $resultCollection =
      CrossmatchSearchEngine::parseOutput( searchOutput => $alignFile );
}
print "Alignment file read in: " . elapsedTime("reading") . "\n";


elapsedTime("scrubbing");
my %consSizeByID = ();
my $numBadRMAlignData = 0;
my $numShort = 0;
my %validate = ();
my %invalid = ();
#open OUT,">tmpSeqList.bed" or die;
my ($tfh, $tfilename) = tempfile("tmpGenSeedsXXXXXXXX", DIR=>".",UNLINK => 0);
for ( my $i = 0 ; $i < $resultCollection->size() ; $i++ ) {
  my $result = $resultCollection->get( $i );
  my $familyName = $result->getSubjName();
  if ( $familyName =~ /\#Simple|\#Low|short/) {
    $invalid{$i} = 1;
    next;
  }
  # Calc cons size
  my $consSize = $result->getSubjEnd() + $result->getSubjRemaining();
  if ( ! exists $consSizeByID{$familyName} ) {
    $consSizeByID{$familyName} = $consSize;
  }elsif ( $consSizeByID{$familyName} != $consSize ) {
    print "WARN: $familyName has more than one reported size: $consSizeByID{$familyName} and $consSize\n";
  }

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
    my $bad = 0;
    $bad = 1 if ( $gapQSeq =~ /.+X/ );

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
  if ( $querySeq =~ /X{10}/ && $options{'assembly'} )
  {
    my $xStart     = $result->getQueryStart() - 1;
    my $twoBitFile = $options{'assembly'};
    my $chr         = $result->getQueryName();
    my $newQuerySeq = "";
    while ( $querySeq =~ /([^X]+)(X{10,}+)/ig )
    {
      my $prefixSeq    = $1;
      my $xSeq         = $2;
      my $prefixSeqLen = $prefixSeq;
      $prefixSeqLen =~ s/-//g;
      $xStart += length( $prefixSeqLen );
      my $xEnd = $xStart + length( $xSeq );
      $newQuerySeq .= $prefixSeq;
      my $replSeq = `$ucscToolsDir/twoBitToFa $twoBitFile:$chr:$xStart-$xEnd stdout`;
      $xStart += length( $xSeq );
      $replSeq =~ s/^>[^\n\r]+[\n\r]+//;
      $replSeq =~ s/[\n\r]//g;
      $newQuerySeq .= $replSeq;
    }
    $newQuerySeq .= substr( $querySeq, length( $newQuerySeq ) );
    $result->setQueryString( $newQuerySeq );
    $querySeq = $newQuerySeq;
  }
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
  if ( length( $qs ) < 30 )
  {
    $numShort++;
    $invalid{$i} = 1;
    next;
  }

  # save twoBitQuery details for sequence validation
  my $twoBitQueryConcat =
         $result->getQueryName() . ":"
      . ( $result->getQueryStart() - 1 ) . "-"
      . $result->getQueryEnd();
  my $twoBitQuery =
         $result->getQueryName() . "\t"
      . ( $result->getQueryStart() - 1 ) . "\t"
      . $result->getQueryEnd() . "\t+\n";
  my $qs = $result->getQueryString();
  $qs =~ s/-//g;
  print $tfh "$twoBitQuery\n";
  $validate{$twoBitQueryConcat} = $qs;
}
close $tfh;

print "$numBadRMAlignData bad RepeatMasker alignment data\n";
print "$numShort sequences were too short ( < 30bp )\n";
print "Scrubbing RM data in: " . elapsedTime("scrubbing") . "\n";
elapsedTime("validating");

#   -bed=input.bed  Grab sequences specified by input.bed. Will exclude introns.
#   -bedPos         With -bed, use chrom:start-end as the fasta ID in output.fa.
open IN,"$ucscToolsDir/twoBitToFa -bedPos -bed=$tfilename $options{'assembly'} stdout|" or die;
my $id = "";
my $seq = "";
my $idx = 0;
my $numFailedSeqValidation = 0;
while ( <IN> ) {
  if ( /^>(\S+)/ ){
    my $tmpId = $1;
    if ( $seq ) {
      if ( exists $validate{$id} ) {
        if ( $validate{$id} ne uc($seq) ) {
          $invalid{$idx} = 1;
          $numFailedSeqValidation++;
          #print "$id did not validate!\n";
        }else {
          # Good mapping
        }
      }else {
        print "WARN: Could not find $id in validate structure\n";
      }
      $idx++;
    }
    $seq = "";
    $id = $tmpId;
    next;
  } 
  s/[\n\r\s]+//g;
  $seq .= $_;
}
if ( $seq ) {
  if ( exists $validate{$id} ) {
    if ( $validate{$id} ne uc($seq) ) {
      $invalid{$idx} = 1;
      $numFailedSeqValidation++;
      #print "$id did not validate!\n";
    }else {
      # Good mapping
    }
  }else {
    print "WARN: Could not find $id in validate structure\n";
  }
}
close IN;
unlink($tfilename);
undef %validate;

my %alignByID = ();
for ( my $i = 0 ; $i < $resultCollection->size() ; $i++ ) {
  unless ( $invalid{$i} ) {
    my $result = $resultCollection->get( $i );
    my $familyName = $result->getSubjName();
    push @{$alignByID{$familyName}}, $result;
  }
}
print "$numFailedSeqValidation sequences failed validation against the assembly.\n";
print "Validating sequences in: " . elapsedTime("validating") . "\n";
undef $resultCollection;
print "  Total Families:  " . scalar( keys( %alignByID ) ) . " ( excluding simple/low )\n";


my $excludeSimple      = 1;
my $excludeShort       = 1;
my $targetSampleCount  = 500;
my $minDepth           = 10;
my $fastaLineLen       = 50;
my $repbaseFile        = "/u1/local/repeatDB/util/reconciliation/data/RepBase20.07/RBCombined.embl";
my $repeatmaskerFile   = "/usr/local/RepeatMasker-4.0.6/Libraries/RepeatMaskerLib.embl";
my $consLen              = 0;
my $totalAttemptedBuilds = 0;
my $assemblyName = "";
# TODO: Make this more general.  
if ( $options{'assembly'} && $options{'assembly'} =~ /(GCA_\d+(\.\d+)?)/ ) {
  $assemblyName = $1;
}
my $DEBUG = 0;


my $totalBuilt = 0;
my $noAlign = 0;
foreach my $id ( keys( %alignByID ) )
{
  $totalAttemptedBuilds++;
  my $countInGenome = scalar(@{$alignByID{$id}});
  $consLen = $consSizeByID{$id};
  print "Working on $id ( length=$consLen, $countInGenome in assembly )\n";
  elapsedTime(1);
  elapsedTime(3);
  
  my @sortedByDiv = sort { $a->getPctKimuraDiverge() <=> $b->getPctKimuraDiverge() } @{$alignByID{$id}};
  my @outliers = ();
  my $divergenceFilter = 100;
  my $medianDiv;
  if ( scalar(@sortedByDiv) % 2 )
  {
    # Odd
    $medianDiv = $sortedByDiv[int(scalar(@sortedByDiv)/2)]->getPctKimuraDiverge();
  }else {
    # Even
    $medianDiv = ($sortedByDiv[int(scalar(@sortedByDiv)/2)-1]->getPctKimuraDiverge() +
               $sortedByDiv[int(scalar(@sortedByDiv)/2)]->getPctKimuraDiverge()
              ) / 2;
  }
  my $quartileDiv = $sortedByDiv[int(scalar(@sortedByDiv)/4)*3]->getPctKimuraDiverge();
  print "  Median Divergence = $medianDiv\n";
  print "  3rd Quartile Divergence = $quartileDiv\n";
  #$divergenceFilter =
  #    ( $medianDiv - $quartileDiv ) + $medianDiv;
  @outliers = splice(@sortedByDiv, int(scalar(@sortedByDiv)/4)*3);
  # Re-sort by length
  @outliers = sort { ($b->getSubjEnd()-$b->getSubjStart()) <=> ($a->getSubjEnd()-$a->getSubjStart()) } @outliers;
  print "  " . scalar(@outliers) . " outliers were moved to the end of the priority list\n";

  # Sort by length, longest first
  my @data = sort { ($b->getSubjEnd()-$b->getSubjStart()) <=> ($a->getSubjEnd()-$a->getSubjStart()) } @sortedByDiv;
  my $outlierStartIdx = scalar(@data);
    
  # Precedence:
  #     Long elements
  #     Elements within first 3 quartiles of kimura divergence
  #     Elements that cover a low-sample-depth region
  push @data, @outliers;

  my $idx         = 0;
  my $sampleCount = 0;
  my $resultCol   = SearchResultCollection->new();
  my @sampledDepth = ();
  my %seen         = ();
  my $numBadReportedConsLen = 0;
  my $numDiagMismatch = 0;
#  print
#"Key: '*' = Saved , '+' = Good Align , '?' = Cons Length , '.' = Bad Align, 'S' = Short, '!' = Incorrect Mapping \n";
  print "   - Sorting and data preparation : " . elapsedTime( 1 ) . "\n";
  elapsedTime(2);
  elapsedTime(4);
  while ( @data )
  {
    my $result = shift @data;
    $idx++;

    # Does this matter?
    my $calcCLen = $result->getSubjEnd() + $result->getSubjRemaining();
    if ( $calcCLen != $consLen )
    {
      #print "?";
      $numBadReportedConsLen++;
    }

    # Must-have positions
    if ( @diagSites )
    {
      my $qry = $result->getQueryString();
      my $sbj = $result->getSubjString();
      
      my @sPosToQBase = ();
      my $sIdx = $result->getSubjStart(); # one based
      $sIdx = $result->setSubjEnd() if ( $result->getOrientation() eq "C" );
      #print "\n$qry\n$sbj\n";
      #print "\n";
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
          #print "  - $sIdx = $qBase\n";
          $sIdx--;  
        }else {
          $sPosToQBase[$sIdx] = $qBase;
          #print "  - $sIdx = $qBase\n";
          $sIdx++;  
        }
      }

      my $failed = 0;
      foreach my $site ( @diagSites )
      {
        if ( $sPosToQBase[$site->[0]] ne $site->[1] )
        {
          #print "  ***FAILED: $site->[0] " . $sPosToQBase[$site->[0]] . " ne " . $site->[1] . "\n"; 
          $failed = 1;
        }
      }

      if ( $failed )
      {
        $numDiagMismatch++;
        #print "&";
        next;
      }
    }

    # Good alignment include it
    #print "+";
    my $qstr = $result->getQueryString();
    $qstr =~ s/X/N/g;
    $result->setQueryString( $qstr );
    $qstr = $result->getSubjString();
    $qstr =~ s/X/N/g;
    $result->setSubjString( $qstr );

    # Remove duplicates! -- hmmm consider this fully
    my $ss  = $result->getSubjStart();
    my $se  = $result->getSubjEnd();
    my $key =
          $result->getScore()
        . $result->getPctDiverge()
        . $result->getPctInsert()
        . $result->getPctDelete();
    next if ( $seen{$key} );
    $seen{$key}++;

    my $addIt           = 0;
    my $covReached      = 1;
    my @samplesToUpdate = ();
    for ( my $i = 0 ; $i < ( $consLen / 10 ) ; $i++ )
    {
      my $idxPos = ( $i * 10 ) + 1;
      $covReached = 0 if ( $sampledDepth[ $i ] < $minDepth );
      next if ( $idxPos > $se );
      next if ( $idxPos < $ss );
      $addIt = 1 if ( $sampledDepth[ $i ] < $minDepth );
      push @samplesToUpdate, $i;
    }
    if ( $sampleCount >= $targetSampleCount && $covReached )
    {
      print "  Coverage reached!\n";
      last;
    }
    if ( ( $idx < $outlierStartIdx && $sampleCount < $targetSampleCount ) || $addIt )
    {
      foreach my $idUpd ( @samplesToUpdate )
      {
        $sampledDepth[ $idUpd ]++;
      }
      #print "*";
      $result->setQueryName( $assemblyName . ":"
                             . $result->getQueryName() ); 
 
      $resultCol->add( $result );
      $sampleCount++;
    }
  } # while ( @data )
  print "   - Selecting elements : " . elapsedTime( 2 ) . "\n";

  my $noCovExamples = 0;
  my $minCovDepth = 10000000000;
  my $maxCovDepth = 0;
  print "  ";
  for ( my $i = 0 ; $i < ( $consLen / 10 ) ; $i++ )
  {
    my $idxPos = ( $i * 10 ) + 1;
    $sampledDepth[$i] = 0 if ( $sampledDepth[$i] eq "" );
    print "[$idxPos]=$sampledDepth[$i], ";
    if ( ($i+1) % 10 == 0 ) {
      print "\n  ";
    }
    $minCovDepth = $sampledDepth[$i] if ( $sampledDepth[$i] < $minCovDepth ); 
    $maxCovDepth = $sampledDepth[$i] if ( $sampledDepth[$i] > $maxCovDepth ); 
    $noCovExamples++ if ( ! defined $sampledDepth[$i]  || $sampledDepth[$i] < 1 );
  }
  print "\n";

  print "  Stats:\n";
  print "    Coverage depth range: $minCovDepth to $maxCovDepth ( from sampled positions )\n";
  print "    $sampleCount of $countInGenome were chosen for the multiple alignment\n";
  print "    $numBadReportedConsLen had differing data about consensus length ( RM artifact )\n";
  print "    $numDiagMismatch had mismatches to the specified diagnostic sites\n";

  if ( $noCovExamples )
  {
    print "  *** Some regions are not covered! ***\n";
    #warn "Some regions of $id are not covered ( $noCovExamples ).\n";
  }

  if ( $minCovDepth < $minDepth )
  {
    print "  *** Some regions did not reach the min coverage depth of $minDepth  ***\n";
    #warn "Some regions of $id are did not reach the min coverage depth of $minDepth.\n";
  }

  if ( $sampleCount == 0 ) {
    print "WARNING: $id is being skipped because there are no alignments?!??\n";
    $noAlign++;
    next;
  }

  my $mAlign = MultAln->new( searchCollection          => $resultCol,
                             searchCollectionReference => MultAln::Subject );

  my $sanitizedID = $id;
  $sanitizedID =~ s/[\(\)]//g; 
  my $class = "Unknown";
  if ( $sanitizedID =~ /(.*)\#(.*)/ )
  {
    $sanitizedID = $1;
    $class = $2;
  }
  
  print "Saving output\n";
  $totalBuilt++;
  if ( $options{'includeRef'} ) {
                      #header           => $newHeaders{uc($id)}->{'header'},
    $mAlign->toSTK(
                      filename         => "$sanitizedID.stk",
                      includeReference => 1,
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
  #my $desc = "Seed alignments generated from RepeatMasker annotations using generateSeedAlignments.pl. ".
  #           "The median Kimura divergence for the family is $medianDiv, $sampleCount where chosen from $countInGenome identified in " . 
  #           "the ". $options{'assembly'} . " assembly file.";
  my $desc = "Source:gsa, mDiv=$medianDiv, $options{'assembly'}:$countInGenome";
  $seedAlign->setDescription($desc);

  if ( $options{'outSTKFile'} ) 
  {
    open OUT,">>".$options{'outSTKFile'} or die;
  }else { 
    open OUT,">$sanitizedID.stk" or die;
  }

  print OUT "" . $seedAlign->toString();
  close OUT;
  print "   - total build time : " . elapsedTime( 3 ) . "\n";
}
print "\n\n";
print "Total Seeds We Attempted To Build: $totalAttemptedBuilds\n";
print "Total built: $totalBuilt\n";
print "Total without alignments: $noAlign\n";
print "\n\n";

# All done
exit;

############################################################################################

##
## Randomize an array
##
sub fisherYatesShuffle
{
  my $array = shift;
  my $i;
  for ( $i = @$array - 1 ; $i >= 0 ; --$i )
  {
    my $j = int rand( $i + 1 );
    next if $i == $j;
    @$array[ $i, $j ] = @$array[ $j, $i ];
  }
}

##-------------------------------------------------------------------------##
##
##  Use: my = writeHTMLMultAlign( multAln => $multAlignRef,
##                                    [destination => $filename|$FH],
##                                    [leftFlankingID => 1],
##                                    [rightFlankingID => 1] );
##
##
##
##-------------------------------------------------------------------------##
my $CLASS = "";

sub writeHTMLMultAlign
{
  my %parameters = @_;

  my $method = "writeHTMLMultAlign";
  croak $CLASS. "::$method() missing multAln parameter!\n"
      if ( !exists $parameters{'multAln'} );

  my $mAlign = $parameters{'multAln'};

  my $OUT = *STDOUT;
  if ( defined $parameters{'destination'} )
  {
    if ( ref( $parameters{'destination'} ) !~ /GLOB|FileHandle/ )
    {
      print $CLASS
          . "::$method() Opening file "
          . $parameters{'destination'} . "\n"
          if ( $DEBUG );
      open $OUT, $parameters{'destination'}
          or die $CLASS
          . "::$method: Unable to open "
          . "results file: $parameters{'destination'} : $!";
    } else
    {
      $OUT = $parameters{'destination'};
    }
  }

  print $OUT <<"END";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
        "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <title>Alignments</title>
  <style type="text/css">
  font.lowQual {
    background: #CD5C5C;
  }
  font.deletion {
    background: #FF0000;
  }
  font.duplication {
    background: #0000FF;
  }
  font.unknown {
    background: #FFFF00;
  }
  font.dupFlank {
    background: #C0C0C0;
  }
  </style>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
</head>
END

  # Find max padding.
  my $maxLeftLen  = 0;
  my $maxRightLen = 0;
  my $maxQueryEnd = length( $mAlign->getReferenceSeq() );
  foreach my $seqNum ( 0 .. $mAlign->getNumAlignedSeqs() - 1 )
  {
    my $relLeftLen =
        length( $mAlign->getLeftFlankingSequence( $seqNum ) ) -
        $mAlign->getAlignedStart( $seqNum );
    my $relRightLen =
        length( $mAlign->getRightFlankingSequence( $seqNum ) ) -
        ( $maxQueryEnd - $mAlign->getAlignedEnd( $seqNum ) );
    $maxLeftLen  = $relLeftLen  if ( $maxLeftLen < $relLeftLen );
    $maxRightLen = $relRightLen if ( $maxRightLen < $relRightLen );
  }

  ## Calculate Del/Dup candidates
  my @dupDelCols = ();
  if ( $parameters{'printDupDelCols'} == 1 )
  {
    my $endStartPairs = $mAlign->_getEndStartPairs();
    foreach my $pair ( @{$endStartPairs} )
    {
      if ( abs( $pair->{'refEnd'} - $pair->{'refStart'} ) < 50 )
      {
        my $adjStart = 0;
        my $adjEnd   = 0;
        my $absStart = $pair->{'refEnd'};
        my $absEnd   = $pair->{'refStart'};
        if ( $pair->{'refEnd'} > $pair->{'refStart'} )
        {
          $absStart = $pair->{'refStart'};
          $absEnd   = $pair->{'refEnd'};
        }
        $adjStart = $mAlign->getAlignPosFromBPPos( $absStart );
        $adjEnd   = $mAlign->getAlignPosFromBPPos( $absEnd );
        print " adjStart = $adjStart adjEnd=$adjEnd\n" if ( $DEBUG );
        if ( $pair->{'avgGapWidth'} > 30 )
        {
          print "Calling a deletion\n";

          # Deletion
          push @dupDelCols,
              {
                'start' => $adjStart,
                'end'   => $adjEnd,
                'type'  => "deletion"
              };

        } elsif ( $pair->{'avgGapWidth'} < 0 )
        {

          # Duplication
          print "Calling a duplication\n";
          my $flnkStart =
              $mAlign->getAlignPosFromBPPos(
                                    $absStart - abs( $pair->{'avgGapWidth'} ) );
          my $flnkEnd = $mAlign->getAlignPosFromBPPos( $absStart - 1 );
          push @dupDelCols,
              {
                'start' => $flnkStart,
                'end'   => $flnkEnd,
                'type'  => 'dupFlank'
              };
          push @dupDelCols,
              {
                'start' => $adjStart,
                'end'   => $adjEnd,
                'type'  => "duplication"
              };
          $flnkStart = $mAlign->getAlignPosFromBPPos( $absEnd + 1 );
          $flnkEnd   =
              $mAlign->getAlignPosFromBPPos(
                                      $absEnd + abs( $pair->{'avgGapWidth'} ) );
          push @dupDelCols,
              {
                'start' => $flnkStart,
                'end'   => $flnkEnd,
                'type'  => 'dupFlank'
              };
        } else
        {

          # Unknown
          print "Calling an Unknown\n";
          push @dupDelCols,
              {
                'start' => $adjStart,
                'end'   => $adjEnd,
                'type'  => "unknown"
              };
        }
      }
    }
  }

  ## Calculate low scoring columns
  my $matrix = SequenceSimilarityMatrix->new();
  $matrix->parseFromFile(
      "$RepModelConfig::REPEATMODELER_MATRICES_DIR/wublast/nt/comparison.matrix"
  );
  my $columns;
  my $scoreArray;
  if ( defined $parameters{'threshold'} )
  {
    ( $columns, $scoreArray ) = $mAlign->getLowScoringAlignmentColumns(
                                           matrix    => $matrix,
                                           threshold => $parameters{'threshold'}
    );
  } else
  {
    ( $columns, $scoreArray ) =
        $mAlign->getLowScoringAlignmentColumns( matrix => $matrix );
  }

  print $OUT "<PRE>\n";

  # Print out scoreArray just for fun
  #   find the largest/smallest score
  my $max = 0;
  for ( my $j = 0 ; $j <= $#{$scoreArray} ; $j++ )
  {
    my $num = sprintf( "%0.1f", $scoreArray->[ $j ] );
    $max = length( $num ) if ( $max < length( $num ) );
    $scoreArray->[ $j ] = $num;
  }
  my $numCols = $max;
  my @lines   = ();
  foreach my $num ( @{$scoreArray} )
  {
    my $paddedNum = ' ' x ( $numCols );
    if ( $num != 0 )
    {
      $paddedNum = ' ' x ( $numCols - length( $num ) ) . $num;
    }
    for ( my $j = 0 ; $j < $numCols ; $j++ )
    {
      $lines[ $j ] .= substr( $paddedNum, $j, 1 );
    }
  }
  my $label = "lowQualScore";
  my $paddedLabel = $label . " " x ( 30 - length( $label ) );
  foreach my $line ( @lines )
  {
    print $OUT $paddedLabel . ": " . ' ' x ( $maxLeftLen ) . $line . "\n";
  }

  # First print the consensus sequence
  my $lineStart = 0;
  my $name      = "consensus";
  my $namePad   = 30 - length( $name );
  my $seq       =
      $mAlign->consensus(
                    "$RepModelConfig::REPEATMODELER_MATRICES_DIR/linupmatrix" );
  print $OUT "<b><i>$name</i></b>"
      . ' ' x $namePad . ": "
      . ' ' x ( $maxLeftLen )
      . "<font color=\"blue\">$seq</font>\n";

  my $name    = "Reference ( " . $mAlign->getReferenceName() . " )";
  my $namePad = 30 - length( $name );
  my $seq     = $mAlign->getReferenceSeq();
  print $OUT "<b><i>$name</i></b>"
      . ' ' x $namePad . ": "
      . ' ' x ( $maxLeftLen )
      . "<font color=\"blue\">$seq</font>\n";

  # Now print the reference and the instances.
  for ( my $i = 0 ; $i < $mAlign->getNumAlignedSeqs() ; $i++ )
  {
    $name = $mAlign->getAlignedName( $i );
    if ( length( $name ) > 30 )
    {
      $name = substr( $name, 0, 30 );
    }
    $namePad = 30 - length( $name );

    my $lfSeq = $mAlign->getLeftFlankingSequence( $i );

    my $seq .=
        ' ' x
        ( $maxLeftLen - ( length( $lfSeq ) - $mAlign->getAlignedStart( $i ) ) );
    if ( $parameters{'leftFlankingID'} eq $i - 1 )
    {
      $seq .= "<font color=\"blue\">" . lc( $lfSeq ) . "</font>";
    } else
    {
      $seq .= lc( $lfSeq );
    }

    # Highlight low scoring columns
    if ( $parameters{'printDupDelCols'} == 1 )
    {
      $seq .= "<b>";
      my $seqPos = 0;
      foreach my $col ( @dupDelCols )
      {
        my $start = $col->{'start'} - $mAlign->getAlignedStart( $i );
        my $end   = $col->{'end'} - $mAlign->getAlignedStart( $i );
        next if ( $start < 0 && $end < 0 );
        $start = 0 if ( $start < 0 );
        $end = length( $mAlign->getAlignedSeq( $i ) ) - $start - 1
            if ( $end < 0 );
        $seq .=
            substr( $mAlign->getAlignedSeq( $i ), $seqPos, $start - $seqPos );
        $seq .= "<font class=\"" . $col->{'type'} . "\">";
        $seq .=
            substr( $mAlign->getAlignedSeq( $i ), $start, $end - $start + 1 );
        $seq .= "</font>";
        $seqPos = $end + 1;
      }
      if ( $seqPos < length( $mAlign->getAlignedSeq( $i ) ) - 1 )
      {
        $seq .= substr( $mAlign->getAlignedSeq( $i ), $seqPos );
      }
      $seq .= "</b>";
    } elsif ( $#{$columns} >= 0 )
    {
      $seq .= "<b>";
      my $seqPos = 0;
      foreach my $col ( @{$columns} )
      {
        my $start = $col->[ 0 ] - $mAlign->getAlignedStart( $i );
        my $end   = $col->[ 1 ] - $mAlign->getAlignedStart( $i );
        next if ( $start < 0 && $end < 0 );
        $start = 0 if ( $start < 0 );
        $end = length( $mAlign->getAlignedSeq( $i ) ) - $start - 1
            if ( $end < 0 );
        $seq .=
            substr( $mAlign->getAlignedSeq( $i ), $seqPos, $start - $seqPos );
        $seq .= "<font class=\"lowQual\">";
        $seq .=
            substr( $mAlign->getAlignedSeq( $i ), $start, $end - $start + 1 );
        $seq .= "</font>";
        $seqPos = $end + 1;
      }
      if ( $seqPos < length( $mAlign->getAlignedSeq( $i ) ) - 1 )
      {
        $seq .= substr( $mAlign->getAlignedSeq( $i ), $seqPos );
      }
      $seq .= "</b>";
    } else
    {
      $seq .= "<b>" . $mAlign->getAlignedSeq( $i ) . "</b>";
    }

    if ( $parameters{'rightFlankingID'} eq $i - 1 )
    {
      $seq .=
            "<font color=\"blue\">"
          . lc( $mAlign->getRightFlankingSequence( $i ) )
          . "</font>";
    } else
    {
      $seq .= lc( $mAlign->getRightFlankingSequence( $i ) );
    }

    print $OUT "<b>$name</b>" . ' ' x $namePad . ": $seq\n";
  }

  if ( defined $parameters{'printHistogram'} )
  {
    print $OUT "\n\n";
    my @columnSeqs = ();
    my $maxRows    = 0;
    foreach my $col ( @{$columns} )
    {
      my ( $cons, $seqsRef ) = $mAlign->getAlignmentBlock(
                                                           start => $col->[ 0 ],
                                                           end   => $col->[ 1 ],
                                                           rawSequences => 1
      );
      push @columnSeqs, [ sort { length( $b ) <=> length( $a ) } @{$seqsRef} ];
      $maxRows = $#{$seqsRef} + 1 if ( $maxRows < ( $#{$seqsRef} + 1 ) );
    }
    $label = "blockSeqs";
    $paddedLabel = $label . " " x ( 30 - length( $label ) );
    for ( my $i = 0 ; $i < $maxRows ; $i++ )
    {
      my $seq = " " x ( $maxLeftLen );
      my $pos = 0;
      for ( my $j = 0 ; $j <= $#{$columns} ; $j++ )
      {
        my $col      = $columns->[ $j ];
        my $colWidth = $col->[ 1 ] - $col->[ 0 ] + 1;
        $seq .= " " x ( $col->[ 0 ] - $pos );
        my $colSeqArray = $columnSeqs[ $j ];
        my $colSeq      = "";
        if ( $#{$colSeqArray} >= 0 )
        {
          $colSeq = shift @{$colSeqArray};
          $colSeq = "." if ( $colSeq eq "" );
        }
        $seq .= $colSeq . " " x ( $colWidth - length( $colSeq ) );
        $pos = $col->[ 1 ] + 1;
      }
      print $OUT "$paddedLabel: $seq\n";
    }
    print $OUT "\n\n";
    if ( defined $parameters{'newConsBlocks'} )
    {
      my $newConsBlocks = $parameters{'newConsBlocks'};
      $label = "blockSeqCons";
      $paddedLabel = $label . " " x ( 30 - length( $label ) );
      my $seq = " " x ( $maxLeftLen );
      my $pos = 0;
      for ( my $j = 0 ; $j <= $#{$newConsBlocks} ; $j++ )
      {
        my $colWidth =
            $newConsBlocks->[ $j ]->{'end'} -
            $newConsBlocks->[ $j ]->{'start'} + 1;
        $seq .= " " x ( $newConsBlocks->[ $j ]->{'start'} - $pos );
        my $colSeq = $newConsBlocks->[ $j ]->{'cons'};
        $seq .=
              "<font color=\"blue\">" . $colSeq
            . "</font>"
            . "*" x ( $colWidth - length( $colSeq ) );
        $pos = $newConsBlocks->[ $j ]->{'end'} + 1;
      }
      print $OUT "$paddedLabel: $seq\n";
    }
    if ( defined $parameters{'finalConsensus'} )
    {
      $label = "originalCons";
      $paddedLabel = $label . " " x ( 30 - length( $label ) );
      my $seq =
            " " x ( $maxLeftLen )
          . "<font color=\"red\">"
          . $mAlign->consensus(
                     "$RepModelConfig::REPEATMODELER_MATRICES_DIR/linupmatrix" )
          . "</font>";
      print $OUT "$paddedLabel: $seq\n";
      $label = "finalCons";
      $paddedLabel = $label . " " x ( 30 - length( $label ) );
      my $seq =
            " " x ( $maxLeftLen )
          . "<font color=\"blue\">"
          . $parameters{'finalConsensus'}
          . "</font>";
      print $OUT "$paddedLabel: $seq\n";
    }
  }

  print $OUT "</PRE>\n";

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
  my %rmToDfamClass = (
    'dna/crypton' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton',
    'dna/crypton-a' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-A',
    'dna/crypton-c' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-C',
    'dna/crypton-f' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-F',
    'dna/crypton-h' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-H',
    'dna/crypton-i' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-I',
    'dna/crypton-r' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-R',
    'dna/crypton-s' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-S',
    'dna/crypton-v' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-V',
    'dna/crypton-x' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Circular_dsDNA_Intermediate;Crypton;Crypton-X',
    'dna/maverick' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;DNA_Polymerase;Maverick',
    'rc' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Rolling_Circle',
    'rc/helitron' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Rolling_Circle;Helitron-1',
    'rc/helitron-2' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Rolling_Circle;Helitron-2',
    'dna' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat',
    'dna/academ-1' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Academ_Group;Academ-1',
    'dna/academ-2' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Academ_Group;Academ-2',
    'dna/academ-h' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Academ_Group;Academ-H',
    'dna/casposons' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;DNA_pol;Casposon',
    'dna/cmc-chapaev' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;CACTA_element;CMC;Chapaev_group;Chapaev',
    'dna/cmc-chapaev-3' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;CACTA_element;CMC;Chapaev_group;Chapaev-3',
    'dna/cmc-enspm' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;CACTA_element;CMC;EnSpm',
    'dna/cmc-mirage' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;CACTA_element;CMC;Mirage',
    'dna/cmc-transib' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;CACTA_element;Transib',
    'dna/dada' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Dada',
    'dna/ginger' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Ginger',
    'dna/hat' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element',
    'dna/hat-ac' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;Activator',
    'dna/hat-blackjack' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;Blackjack',
    'dna/hat-charlie' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;Charlie',
    'dna/hat-pegasus' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;Pegasus',
    'dna/hat-restless' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;Restless',
    'dna/hat-tag1' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;Tag1',
    'dna/hat-tip100' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;Tip100',
    'dna/hat-hat1' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hAT1',
    'dna/hat-hat19' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hAT19',
    'dna/hat-hat5' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hAT5',
    'dna/hat-hat6' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hAT6',
    'dna/hat-hatm' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hATm',
    'dna/hat-hatw' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hATw',
    'dna/hat-hatx' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hATx',
    'dna/hat-hobo' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;hAT_Element;hobo',
    'dna/is3eu' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;IS3EU',
    'dna/kolobok' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Kolobok_Group;Kolobok',
    'dna/kolobok-e' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Kolobok_Group;Kolobok-E',
    'dna/kolobok-h' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Kolobok_Group;Kolobok-H',
    'dna/kolobok-hydra' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Kolobok_Group;Hydra-specific_Branch',
    'dna/kolobok-t2' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Kolobok_Group;T2',
    'dna/mule' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Mutator-like_Element',
    'dna/mule-f' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Mutator-like_Element;F',
    'dna/mule-mudr' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Mutator-like_Element;MuDR',
    'dna/mule-nof' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Mutator-like_Element;NOF',
    'dna/merlin' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Merlin',
    'dna/novosib' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Novosib',
    'dna/p' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;P_Element',
    'dna/p-fungi' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;P_Element;Fungi-specific_Branch',
    'dna/pif-harbs' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;PIF-like_Elements;HarbS',
    'dna/pif-harbinger' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;PIF-like_Elements;Harbinger',
    'dna/pif-isl2eu' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;PIF-like_Elements;ISL2EU',
    'dna/pif-spy' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;PIF-like_Elements;Spy',
    'dna/piggybac' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;PiggyBac-like_element;PiggyBac',
    'dna/piggybac-a' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;PiggyBac-like_element;PiggyBac-A',
    'dna/piggybac-x' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;PiggyBac-like_element;PiggyBac-X',
    'dna/sola-1' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Sola-group;Sola-1',
    'dna/sola-2' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Sola-group;Sola-2',
    'dna/sola-3' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Sola-group;Sola-3',
    'dna/tcmar-ant1' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Ant1',
    'dna/tcmar-cweed' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Cweed',
    'dna/tcmar-gizmo' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Gizmo',
    'dna/tcmar-isrm11' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;ISRm11',
    'dna/tcmar-m44' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;m44',
    'dna/tcmar-mariner' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Mariner',
    'dna/tcmar-mogwai' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Mogwai',
    'dna/tcmar-fot1' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Pogo-group;Fot1',
    'dna/tcmar-pogo' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Pogo-group;Pogo',
    'dna/tcmar-tigger' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Pogo-group;Tigger',
    'dna/tcmar-sagan' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Sagan',
    'dna/tcmar-stowaway' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Stowaway',
    'dna/tcmar-tc1' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Tc1',
    'dna/tcmar-tc2' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Tc2',
    'dna/tcmar-tc4' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element;Tc4',
    'dna/tcmar' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Tc1-Mariner-like_Element',
    'dna/zator' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Zator',
    'dna/zisupton' =>
'Interspersed_Repeat;Transposable_Element;DNA_Transposon;Terminal_Inverted_Repeat;Zisupton',
    'line' =>
        'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE',
    'line/cre-ambal' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-I;Ambal',
    'line/cre' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-I;CRE
',
    'line/cre-1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-I;CRE;CRE-1',
    'line/cre-2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-I;CRE;CRE-2',
    'line/cre-odin' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-I;Odin',
    'line/genie' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Genie',
    'line/l1-dre' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;L1-like;DRE',
    'line/l1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;L1-like;L1-group;L1',
    'line/l1-tx1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;L1-like;L1-group;Tx1',
    'line/l1-zorro' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;L1-like;Zorro',
    'line/proto1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;Proto-1',
    'line/r2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;R2-like;R2',
    'line/r2-hero' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;R2-like;Hero',
    'line/r2-nesl' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;R2-like;NeSL',
    'line/deceiver' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;R4-like;Deceiver',
    'line/dong-r4' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;R4-like;Dong-R4',
    'line/dualen' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-1;R4-like;Dualen',
    'line/proto2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;Proto-2',
    'line/cr1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;CR1-group;CR1',
    'line/cr1-zenon' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;CR1-group;CR1;Zenon',
    'line/l2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;CR1-group;L2-group;L2',
    'line/rex-babar' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;CR1-group;Rex-Babar',
    'line/i' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;R1-group;I-group;I',
    'line/i-jockey' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;R1-group;I-group;Jockey',
    'line/r1-loa' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;R1-group;R1-subgroup;LOA',
    'line/r1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;R1-group;R1-subgroup;R1',
    'line/tad1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;R1-like;Tad1',
    'line/rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;RTE-like',
    'line/rte-bovb' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;RTE-like;RTE-group;BovB',
    'line/rte-orte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;RTE-like;ORTE',
    'line/rte-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;RTE-like;RTE-group;RTE',
    'line/rte-x' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE;Group-II;Group-2;RTE-like;RTE-group;RTE-X',
    'retroposon' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter',
    'sine/i' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;I-derived',
    'retroposon/sva' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-dependent;SVA',
    'retroposon/l1-dep' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-dependent',
    'retroposon/rte-derived' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;RTE-derived',
    'sine/l1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-derived',
    'sine/l2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L2-derived',
    'sine/dong-r4' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;R4-derived',
    'sine' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE',
    'sine/5s-deu' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Deu-core;Unknown_LINE-dependent',
    'sine/5s-deu-l2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Deu-core;L2-end',
    'sine/5s-core-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;MIR-core;RTE-end',
    'sine/5s-sauria-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Sauria-core;RTE-end',
    'sine/5s' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter',
    'sine/5s-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;No_or_Unknown_Core;RTE-end',
    'sine/alu' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter;No-core;L1-dependent;Alu',
    'sine/b2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter;No-core;L1-dependent;B2',
    'sine/7sl' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter',
    'sine/trna-5s' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_and_5S_RNA;No_or_Unknown_Core;Unknown_LINE-dependent',
    'sine/b4' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_and_7SL_RNA;No-core;L1-dependent',
    'sine/trna-7sl' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_and_7SL_RNA;No_or_Unknown_Core;Unknown_LINE-dependent',
    'sine/trna-ceph' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Ceph-core;Unknown_LINE-dependent',
    'sine/trna-ceph-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Ceph-core;RTE-end',
    'sine/trna-deu' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;Unknown_LINE-dependent',
    'sine/trna-deu-cr1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;CR1-end',
    'sine/trna-deu-i' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;I-end',
    'sine/trna-deu-l2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;L2-end',
    'sine/trna-deu-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;RTE-end',
    'sine/trna-meta' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Meta-core;Unknown_LINE-dependent',
    'sine/mir' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;L2-end',
    'sine/trna-core' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;Unknown_LINE-dependent',
    'sine/trna-core-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;RTE-end',
    'sine/trna-mermaid' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;Mermaid',
    'sine/id' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No-core;L1-dependent',
    'sine/rte-bovb' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;BovB-end',
    'sine/trna-cr1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;CR1-end',
    'sine/trna-i' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;I-end',
    'sine/trna-jockey' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Jockey-end',
    'sine/trna-l1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;L1-dependent',
    'sine/trna-l2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;L2-end',
    'sine/r1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;R1-end',
    'sine/trna-r2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;R2-end',
    'sine/trna-rex' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Rex-end',
    'sine/trna-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;RTE-end',
    'sine/trna-tad1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Tad1_End',
    'sine/trna-sauria' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;Unknown_LINE-dependent',
    'sine/trna-sauria-l2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;L2-end',
    'sine/trna-sauria-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;RTE-end',
    'sine/trna-v-core-l2' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V_and_MIR-core;L2-end',
    'sine/trna-v' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;Unknown_LINE-dependent',
    'sine/trna-v-cr1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;CR1-end',
    'sine/trna' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Unknown_LINE-dependent',
    'sine/u' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;U-RNA_Promoter;No_or_Unknown_Core;Unknown_LINE-dependent',
    'sine/ceph' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;Unknown_Promoter;Ceph-core;RTE-end',
    'sine/core' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;Unknown_Promoter;MIR-core;Unknown_LINE-dependent',
    'sine/core-rte' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;LINE-dependent_Retroposon;SINE;Unknown_Promoter;MIR-core;RTE-end',
    'ltr/dirs' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Inverted_Long_Terminal_Repeat_Elements;Tyrosine_Recombinase_Elements;DIRS',
    'ltr/ngaro' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Inverted_Long_Terminal_Repeat_Elements;Tyrosine_Recombinase_Elements;Ngaro',
    'ltr/viper' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Inverted_Long_Terminal_Repeat_Elements;Tyrosine_Recombinase_Elements;Viper',
    'ltr' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element',
    'ltr/pao' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Bel-Pao',
    'ltr/gypsy' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Gypsy',
    'ltr/caulimovirus' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Pararetroviridae;Caulimoviridae',
    'ltr/erv1' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1',
    'ltr/erv-lenti' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2;Lenti',
    'ltr/ervk' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV2',
    'ltr/ervl' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV3;ERVL
',
    'ltr/ervl-malr' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV3;MaLR
',
    'ltr/erv4' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV4',
    'ltr/erv' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae',
    'ltr/erv-foamy' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Spumaretrovirinae',
    'ltr/cassandra' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae',
    'ltr/trim' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;TRIM',
    'ltr/copia' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Long_Terminal_Repeat_Element;Ty1-Copia',
    'line/penelope' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;Penelope-like_Elements',
    'unknown/tate' =>
'Interspersed_Repeat;Transposable_Element;Retrotransposed_Element;Retrotransposon;TATE',
    'rrna'                   => 'Interspersed_Repeat;Pseudogene;RNA;rRNA',
    'scrna'                  => 'Interspersed_Repeat;Pseudogene;RNA;scRNA',
    'snrna'                  => 'Interspersed_Repeat;Pseudogene;RNA;snRNA',
    'trna'                   => 'Interspersed_Repeat;Pseudogene;RNA;tRNA',
    'unknown'                => 'Interspersed_Repeat;Unknown',
    'unknown/centromeric'    => 'Interspersed_Repeat;Unknown;Centromeric',
    'satellite'              => 'Tandem_Repeat;Satellite',
    'satellite/acromeric'    => 'Tandem_Repeat;Satellite;Acromeric',
    'satellite/centromeric'  => 'Tandem_Repeat;Satellite;Centromeric',
    'satellite/macro'        => 'Tandem_Repeat;Satellite;Macro',
    'satellite/subtelomeric' => 'Tandem_Repeat;Satellite;Subtelomeric',
    'satellite/w-chromosome' => 'Tandem_Repeat;Satellite;W-chromosomal',
    'satellite/y-chromosome' => 'Tandem_Repeat;Satellite;Y-chromosomal',
    'simple_repeat'          => 'Tandem_Repeat;Simple',
    'other/dna_virus'        => 'Accidental;Normally_Non-integrating_Virus',
    'artefact'               => 'Artifact ',
    'low_complexity'         => 'Low_Complexity',
    'other'                  => 'Other',
    'segmental'              => 'Segmental_Duplication',
  );

  if ( exists $rmToDfamClass{$rmClass} ){
    return $rmToDfamClass{$rmClass};
  }else {
    return "Unknown";
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
    my $DiffTime = time - $TimeBefore{$TimeHistIdx};
    $TimeBefore{$TimeHistIdx} = time;
    my $Min = int( $DiffTime / 60 );
    $DiffTime -= $Min * 60;
    my $Hours = int( $Min / 60 );
    $Min -= $Hours * 60;
    my $Sec = $DiffTime;
    my $timeStr = sprintf( "%02d:%02d:%02d", $Hours, $Min, $Sec );
    return "$timeStr (hh:mm:ss) Elapsed Time";
  }
  else {
    $TimeBefore{$TimeHistIdx} = time;
    return 0;
  }
}

1;
