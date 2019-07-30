#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) viewMSA
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility program to assist with viewing a small multiple alignment
##      using HTML.
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
# ChangeLog
#
#     $Log: viewMSA.pl,v $
#     Revision 1.16  2017/04/05 00:03:32  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

viewMSA - View small multiple alignments using HTML

=head1 SYNOPSIS

  viewMSA [-version] [-order <start|div>] [-useOld]
                      -sequences <*.fa> -ref_sequence <*.fa>
                      [-clustalTree <*.phb> ]

or

  viewMSA [-version] [-order <start|div>] [-useOld]
                      -search_results <cross_match>
                      [-clustalTree <*.phb> ]

or

  viewMSA [-version] [-useOld]
                      -malign *.malign



=head1 DESCRIPTION

  A utility script to view small multiple alignments using HTML.  Strictly
speaking this is really a one-vs-all alignment and not a comprehensive
multiple alignment. 

The options are:

=over 4

=item -sequences <*.fa>

The file containing all the sequences to align to the reference.

=item -ref_sequence <*.fa>

The file containing the reference sequence.

=item -search_results <crossmatch file>

View a one-vs-all alignment file in cross_match format. 

=item -order <start|div>

Order alignment instances by start position or by Kimura
divergence.

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2012-2014 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use JSON::PP;
use Carp;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/..";
use RepModelConfig;

# RepeatMasker Libraries
use lib $RepModelConfig::REPEATMASKER_DIR;
use MultAln;
use SequenceSimilarityMatrix;
use CrossmatchSearchEngine;
use FastaDB;
use SeqDBI;

#
# Hopefully crossmatch is defined here
use RepeatMaskerConfig;

#
# Version
#  -- NOTE: This is filled in by configure
my $Version = "open-1.0.11";
$Version = "DEV" if ( $Version =~ /\#VERSION\#/ );

##----------------------------------------------------------------------##
##      S I T E   S P E C I F I C   C O N F I G U R A T I O N
##
##  If you must include site specific variables in the program
##  itself put them here.
##
##  ie. my $blastPrgrmDir = "/user/local/blast/bin";
##      my $indelPenalty = 30;
##
##  END OF SITE SPECIFIC CONFIGURATION -- DO NOT EDIT BELOW THIS LINE
##----------------------------------------------------------------------##

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',            # print out the version and exit
                    '-sequences=s',
                    '-ref_sequence=s',
                    '-search_results=s',
                    '-useOld',
                    '-malign=s',
                    '-order=s',
                    '-clustalTree=s',
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) )
{
  usage();
}

sub usage
{
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} )
{
  print "$Version\n";
  exit;
}

my $searchResultsFile;
my $order;
if ( $options{'order'} )
{
  if ( $options{'order'} =~ /start|div/i )
  {
    $order = lc( $options{'order'} );
  } else
  {
    print "\n\nOrder parameter not recognized ( $options{'order'} ). "
        . "Must be either 'start' or 'div'!\n\n";
    usage();
  }
}

my $resultCollection;
my $refInput = MultAln::Subject;

if ( $options{'search_results'} )
{
  $searchResultsFile = $options{'search_results'};
  $resultCollection  =
      CrossmatchSearchEngine::parseOutput( searchOutput => $searchResultsFile );

  # TODO: Deprecate this and move it to SearchResultCollection.pm
  # Auto detect which input ( query/subject ) is the static sequence for
  # which all other sequences are aligned.
  my $queryID;
  my $subjID;
  my $staticQuery   = 1;
  my $staticSubject = 1;
  for ( my $i = 0 ; $i < $resultCollection->size() ; $i++ )
  {
    my $result = $resultCollection->get( $i );
    my $qID    = $result->getQueryName();
    my $sID    = $result->getSubjName();
    $staticQuery   = 0 if ( defined $queryID && $queryID ne $qID );
    $staticSubject = 0 if ( defined $subjID  && $subjID  ne $sID );
    die "Strange...this appears not to be a multiple alignment!"
        if ( $staticQuery == 0 && $staticSubject == 0 );
    $queryID = $qID;
    $subjID  = $sID;
  }
  die "Could not determine reference sequence.  This doesn't look like\n"
      . "a multiple alignment to one reference sequence!\n"
      if ( $staticQuery && $staticSubject );
  $refInput = MultAln::Query if ( $staticQuery );

} elsif ( !$options{'malign'} )
{
  my $refSeqFile;
  my $compSeqFile;
  if ( !$options{'sequences'} )
  {
    print "\nMissing option -sequences <*.fa>!\n\n";
    usage();
  }
  $compSeqFile = $options{'sequences'};
  if ( !$options{'ref_sequence'} )
  {
    print "\nMissing option -ref_sequence <*.fa>!\n\n";
    usage();
  }
  $refSeqFile = $options{'ref_sequence'};

  if ( !-s $RepeatMaskerConfig::CROSSMATCH_PRGM )
  {
    die "Could not find the crossmatch program.  Perhaps RepeatModeler\n"
        . "hasn't been configured yet?\n";
  }
  my $cmCmd =
        "$RepeatMaskerConfig::CROSSMATCH_PRGM $compSeqFile $refSeqFile "
      . "-matrix $RepModelConfig::REPEATMODELER_DIR/Matrices"
      . "/crossmatch/comparison.matrix "
      . "-gap_init -25 -del_gap_ext -5 -ins_gap_ext -5 "
      . "-minscore 150 -minmatch 7 -alignments -masklevel 80  2> /dev/null "
      . "> tmpResults.out ";
  print "Running comparison...\n";
  `$cmCmd`;
  $resultCollection =
      CrossmatchSearchEngine::parseOutput( searchOutput => "tmpResults.out" );

  ## TODO: Add some error checking here...

  #unlink( "tmpResults.out" )      if ( -e "tmpResults.out" );
  unlink( "refined_subs.fa.log" ) if ( -e "refined_subs.fa.log" );
}

my $malign;
if ( $resultCollection && $resultCollection->size() > 0 )
{
  print "  -- Building multiple alignment\n";
  print "Search result collection = " . $resultCollection->size() . "\n";

  if ( $options{'clustalTree'} )
  {
    my @idOrder = loadClustalWTreeProximity( $options{'clustalTree'} );
    my %idOrderHash;
    for ( my $i = 0 ; $i <= $#idOrder ; $i++ )
    {
      $idOrderHash{ $idOrder[ $i ] } = $i;
    }
    $resultCollection->sort(
      sub ($$) {
        ( $idOrderHash{ $_[ 0 ]->getQueryName() } <=>
          $idOrderHash{ $_[ 1 ]->getQueryName() } );
      }
    );
  }
  $malign = MultAln->new( searchCollection          => $resultCollection,
                          searchCollectionReference => $refInput, );
} elsif ( $options{'malign'} )
{
  $malign = MultAln->new();
  $malign = $malign->serializeIN( $options{'malign'} );
} else
{
  die "Could not build multiple alignment -- no result collection!\n";
}

if ( $malign )
{
  my $MOUT;
  open $MOUT, ">MultipleAlignment.html";

  if ( $options{'useOld'} )
  {
    print $MOUT "<HR><P><H1>Heading...</H1>\n";
    &writeHTMLMultAlign(
                         multAln              => $malign,
                         destination          => $MOUT,
                         noteTransitionsWithI => 1,
                         highlightCpGs        => 1
    );
  } else
  {
    generateJavascriptSummaryAndAlignmentViewer( $malign, $MOUT );
  }

  close $MOUT;
}

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
##
##  Use: my = writeHTMLMultAlign( multAln => $multAlignRef,
##                                [destination => $filename|$FH],
##                                [leftFlankingID => 1],
##                                [rightFlankingID => 1] );
##
##
##
##-------------------------------------------------------------------------##
sub writeHTMLMultAlign
{
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  croak "$subroutine(" . ") missing multAln parameter!\n"
      if ( !exists $parameters{'multAln'} );

  my $mAlign = $parameters{'multAln'};

  my $OUT = *STDOUT;
  if ( defined $parameters{'destination'} )
  {
    if ( ref( $parameters{'destination'} ) !~ /GLOB|FileHandle/ )
    {
      print "::$subroutine("
          . ") Opening file "
          . $parameters{'destination'} . "\n"
          if ( $DEBUG );
      open $OUT, $parameters{'destination'}
          or die
          . "$subroutine("
          . ": Unable to open "
          . "results file: $parameters{'destination'} : $!";
    } else
    {
      $OUT = $parameters{'destination'};
    }
  }

  ## Print the HTML header
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
  font.CpG {
    background: #FAE2E8;
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
END
  print $OUT '<script type="text/javascript">';
  print $OUT inlineTableJS() . "\n";
  print $OUT '</script>';
  print $OUT "</head>\n";
  print $OUT
'<H1>Summary View</H1><button onClick="visualizeMultiple(\'orient\');">Orientation Sort</button>
<button onClick="visualizeMultiple(\'norm\');">Normal Sort</button>
<button onClick="visualizeMultiple(\'end\');">End Sort</button>
<button onClick="visualizeMultiple(\'div\');">Divergence Sort</button>
<br>
<p>
    <div id="canvasesdiv" style="position:relative">
        <canvas id="alignment_canvas" width="800" height="1600" style="z-index:1;position:absolute;left:0px;top:0px;">Canvas not supported</canvas>
        <canvas id="detail_canvas" width="1200" height="2000" style="z-index:2;position:absolute;left:0px;top:0px;">Canvas not supported</canvas>
        <canvas id="guideline_canvas" width="800" height="1600" style="z-index:3;position:absolute;left:0px;top:0px;">Canvas not supported</canvas>
    </div>
';
  print $OUT
"<H1>Detail View</H1><PRE><table id=\"table-2\" align=\"left\" cellspacing=\"0\" cellpadding=\"2\">\n";

  my $referenceNamePrefix = "Reference ( ";
  my $referenceNameSuffix = " )";

  # Get ordered set of seqIDs from the seqDB
  my @orderedSeqIDs  = ();
  my @orderedSeqNums = ();
  my %desc           = ();
  if ( defined $parameters{'SeqDB'} )
  {
    @orderedSeqIDs = $parameters{'SeqDB'}->getIDs();
  }

  # Find max padding & max name/desc length
  my $maxLeftLen        = 0;
  my $maxRightLen       = 0;
  my $maxQueryEnd       = length( $mAlign->getReferenceSeq() );
  my $maxNameDescLength = 0;
  my %nameToNum         = ();
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

    my $name   = $mAlign->getAlignedName( $seqNum );
    my $curLen = length( $name );
    $curLen += length( $referenceNamePrefix ) + length( $referenceNameSuffix )
        if ( $seqNum == 0 );
    $nameToNum{$name} = $seqNum;
    if (    defined $parameters{'useSeqDBDesc'}
         && defined $parameters{'SeqDB'} )
    {
      $desc{$name} = $parameters{'SeqDB'}->getDescription( $name );
      $curLen += length( $desc{$name} ) + 1;
    }
    $maxNameDescLength = $curLen if ( $curLen > $maxNameDescLength );
  }

  # If provided a SeqDB then order the instances by
  # SeqDB
  if ( defined $parameters{'SeqDB'} )
  {
    for ( my $i = 0 ; $i <= $#orderedSeqIDs ; $i++ )
    {
      if ( exists $nameToNum{ $orderedSeqIDs[ $i ] } )
      {
        $orderedSeqIDs[ $i ] = $nameToNum{ $orderedSeqIDs[ $i ] };
      } else
      {
        warn "Could not determine the id number for $orderedSeqIDs[$i]\n";
        $orderedSeqIDs[ $i ] = -1;
      }
    }
    unshift @orderedSeqIDs, 0;
  } else
  {
    @orderedSeqIDs = ( 0 .. $mAlign->getNumAlignedSeqs() );
    if ( $order eq "div" )
    {

      # Order by divergence
      $mAlign->kimuraDivergence();
      @orderedSeqIDs =
          sort { $mAlign->getAlignedDiv( $a ) <=> $mAlign->getAlignedDiv( $b ) }
          ( 0 .. ( $mAlign->getNumAlignedSeqs() - 1 ) );
      unshift( @orderedSeqIDs, 0 );
    } elsif ( $order eq "start" )
    {
      @orderedSeqIDs = sort {
        $mAlign->getAlignedStart( $a ) <=> $mAlign->getAlignedStart( $b )
      } ( 0 .. ( $mAlign->getNumAlignedSeqs() - 1 ) );
      unshift( @orderedSeqIDs, 0 );
    }
  }

  # First print the consensus sequence
  my $lineStart = 0;
  my $name      = "consensus";
  my $namePad   = $maxNameDescLength - length( $name );

  #print "About to call for a consensus\n";
  my $seq = $mAlign->consensus();
  my $ref = $mAlign->getReferenceSeq();

  # Calculate the ungapped reference length
  my $ugRef = $ref;
  $ugRef =~ s/-//g;
  my $ugRefLen = length( $ugRef );
  $ugRef = undef;

  # Print position index
  my @lines = ();
  for ( my $i = 0 ; $i < length( $ugRefLen ) ; $i++ )
  {
    $lines[ $i ] = ' ' x $maxNameDescLength . "   " . ' ' x ( $maxLeftLen );
  }
  my $idxInt = 50;
  my $pos    = 0;
  for ( my $j = 0 ; $j < length( $ref ) ; $j++ )
  {
    my $refChar = substr( $ref, $j, 1 );
    my $idxLen = 0;
    if ( $refChar ne "-" )
    {
      $pos++;
      if ( $pos == 1 || $pos % $idxInt == 0 || $pos == $ugRefLen )
      {
        $idxLen = length( $pos );
      }
    }
    for ( my $k = 0 ; $k < length( $ugRefLen ) ; $k++ )
    {
      if ( $k < $idxLen )
      {
        $lines[ $k ] .= substr( $pos, length( $pos ) - $k - 1, 1 );
      } else
      {
        $lines[ $k ] .= " ";
      }
    }
  }
  print $OUT "<pre>\n";
  for ( my $i = $#lines ; $i >= 0 ; $i-- )
  {
    print $OUT $lines[ $i ] . "\n";
  }
  print $OUT "</pre>\n";

  # Printing Consensus
  print $OUT
      "<tr><th align=\"left\" nowrap=\"nowrap\"><pre><b><i>$name</i></b> "
      . ' ' x $namePad . ": "
      . ' ' x ( $maxLeftLen );

  warn "Ref/cons not the same length!\n"
      if ( length( $ref ) != length( $seq ) );
  for ( my $u = 0 ; $u < length( $ref ) ; $u++ )
  {
    if ( substr( $seq, $u, 1 ) ne substr( $ref, $u, 1 ) )
    {
      print $OUT "<font color=\"blue\">" . substr( $seq, $u, 1 ) . "</font>";
    } else
    {
      print $OUT "" . substr( $seq, $u, 1 );
    }
  }
  print $OUT "</pre></th></tr>\n";

  # Now print the reference and the instances.
  $name    = "Reference ( " . $mAlign->getReferenceName() . " )";
  $namePad = $maxNameDescLength - length( $name );
  print $OUT "<tr><th align=\"left\" nowrap=\"nowrap\"><pre><b>$name</b> "
      . ' ' x $namePad
      . ": $ref</pre></th></tr>\n";

  foreach my $i ( @orderedSeqIDs )
  {
    next if $i < 0;
    $name = $mAlign->getAlignedName( $i );
    my $desc = "";
    if (    defined $parameters{'useSeqDBDesc'}
         && defined $parameters{'SeqDB'} )
    {
      $desc = " " . $desc{ $mAlign->getAlignedName( $i ) };
    }
    $namePad = $maxNameDescLength - length( $name . $desc );

    my $lfSeq = $mAlign->getLeftFlankingSequence( $i );

    my $leftPad = 0;
    my $start   = $mAlign->getAlignedStart( $i );
    $leftPad = ( $maxLeftLen - ( length( $lfSeq ) - $start ) );
    $seq = ' ' x ( $leftPad );
    if ( $parameters{'leftFlankingID'} eq $i - 1 )
    {
      $seq .= "<font color=\"blue\">" . lc( $lfSeq ) . "</font>";
    } else
    {
      $seq .= lc( $lfSeq );
    }
    for ( my $j = $mAlign->getAlignedStart( $i ) ;
          $j <= $mAlign->getAlignedEnd( $i ) ;
          $j++ )
    {
      my $base = substr( $mAlign->getAlignedSeq( $i ),
                         $j - $mAlign->getAlignedStart( $i ), 1 );
      my $cb          = substr( $mAlign->getAlignedSeq( 0 ), $j,     1 );
      my $suffixDiNuc = substr( $mAlign->getAlignedSeq( 0 ), $j,     2 );
      my $prefixDiNuc = substr( $mAlign->getAlignedSeq( 0 ), $j - 1, 2 );
      if (
        $parameters{'highlightCpGs'}
        &&

        #$suffixDiNuc =~ /(CG|TG|CA)/ )
        $suffixDiNuc =~ /CG/
          )
      {
        $seq .= "<font class=\"CpG\">";
      }
      if ( $base ne "-" && $base eq $cb && $i > 0 )
      {
        $seq .= ".";
      } elsif ( $base eq "-" && $base eq $cb && $i > 0 )
      {
        $seq .= " ";
      } else
      {
        ##
        ##  Note transitions with an 'i' if desired
        ##    A->G or C->T Transition
        ##
        if (
                $base ne $cb
             && $parameters{'noteTransitionsWithI'}
             && (    $cb =~ /[CT]/ && $base =~ /[CT]/
                  || $cb =~ /[AG]/ && $base =~ /[AG]/ )
            )
        {
          $seq .= "i";
        } else
        {
          $seq .= $base;
        }
      }

      if (
        $parameters{'highlightCpGs'}
        &&

        #$prefixDiNuc =~ /(CG|TG|CA)/ )
        $prefixDiNuc =~ /CG/
          )
      {
        $seq .= "</font>";
      }

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

    print $OUT
        "<tr><th align=\"left\" nowrap=\"nowrap\"><pre><b>$name</b>$desc "
        . ' ' x $namePad
        . ": $seq</pre></th></tr>\n";
  }

  print $OUT "</table></PRE>\n";
  print $OUT "<script type=\"text/javascript\">\n";
  print $OUT "var table = document.getElementById('table-2');\n";
  print $OUT "var tableDnD = new TableDnD();\n";
  print $OUT "tableDnD.init(table);\n";
  print $OUT "</script>\n";
  print $OUT '<script type="text/javascript">';
  print $OUT "\nvar data = " . &getMSASummaryViewJSON( $mAlign ) . ";\n";
  print $OUT summaryViewJS() . "\n";
  print $OUT '</script>';
}

#
# Load in a ClustalW *.phb file ( tree ) and use
# proximity of ids in that file as a way to order
# the ids in the multiple alignment.
#
sub loadClustalWTreeProximity
{
  my $fileName = shift;

  open IN, "<$fileName"
      or die "Could not open Clustalw tree file ( $fileName ): $!\n";
  my @seqIDOrder = ();
  while ( <IN> )
  {

    # (
    # (
    # name1:distance[bootstrap], name2:distance[bootstrap])
    # :distance[bootstrap]
    # );
    # bootstrap is not always present
    if ( /(\S+)\:[\d\.]+/ )
    {
      my $id = $1;
      $id =~ s/\(\)//g;
      push @seqIDOrder, $id;
    }
  }
  close IN;
  return ( @seqIDOrder );
}

sub summaryViewJS
{
  my $javascript = '
//
// Setup HTML5 Canvas
//
var align_canvas = document.getElementById("alignment_canvas");
var guide_canvas = document.getElementById("guideline_canvas");
align_canvas.height = ( data.num_alignments * 2 ) + 10 + 8 + 10 + 20;
guide_canvas.height = align_canvas.height;
var cdiv = document.getElementById("canvasesdiv")
cdiv.style.height = align_canvas.height;
// Get drawing contexts
var guide_context = guide_canvas.getContext("2d");
var align_context = align_canvas.getContext("2d");

// Heatmap Colors
var qualColor = ["#ff6600", "#ffcc00", "#ccff00", "#66ff00", "#00ff00",
    "#00ff66", "#00ffcc", "#00ccff", "#0066ff", "#0000ff"];

// Constants to reduce lookup(?) in event listener
var WIDTH = align_canvas.width;
var HEIGHT = align_canvas.height;
var pixelToBP = data.length / WIDTH;
var currRulerY = 0;


function getMousePos(canvas, evt) {
    var rect = canvas.getBoundingClientRect();
    return {
        x: evt.clientX - rect.left,
        y: evt.clientY - rect.top
    };
}


guide_canvas.addEventListener("mousemove", function (evt) {
    var mousePos = getMousePos(guide_canvas, evt);
    if (mousePos.x >= 10) {
        guide_context.clearRect(0, 0, WIDTH, HEIGHT);
        guide_context.strokeStyle = "#ff0000";
        guide_context.beginPath();
        guide_context.moveTo(mousePos.x, 0);
        guide_context.lineTo(mousePos.x, HEIGHT);
        guide_context.stroke();

        guide_context.font = "italic 11pt Calibri";
        var txt = "" + Math.round(((mousePos.x - 10) * pixelToBP) + 1);
        var text_width = guide_context.measureText(txt).width;
        var text_height = 12; //Estimated based on font ( no height call in HTML5 )
        var textXPos = mousePos.x - (text_width / 2);
        // TODO: Use coordinates of ruler....get somehow
        if (textXPos < 10) {
            textXPos = 10;
        }
        if (textXPos + text_width > WIDTH) {
            textXPos = WIDTH - text_width - 10;
        }

        guide_context.fillStyle = "#FAF7F8";
        guide_context.fillRect(textXPos, currRulerY, text_width, text_height);
        guide_context.fillStyle = "#000000";

        guide_context.fillText(txt, textXPos, currRulerY + 11);

    }
}, false);


//
// 
//
function ruler(x, y, width, height, minVal, maxVal, minorTickInterval, majorTickInterval) {
    align_context.beginPath();
    align_context.moveTo(x, y);
    align_context.lineTo(x + width, y);
    align_context.stroke();

    // Translations
    var pixelsPerUnit = width / (maxVal - minVal + 1);
    var pixelsPerMajorTick = majorTickInterval * pixelsPerUnit;
    var pixelsPerMinorTick = minorTickInterval * pixelsPerUnit;

    for (var i = 0; i < width; i += pixelsPerMajorTick) {
        align_context.beginPath();
        align_context.moveTo(x + i, y);
        align_context.lineTo(x + i, y + height);
        align_context.stroke();
    }

    for (var i = 0; i < width; i += pixelsPerMinorTick) {
        align_context.beginPath();
        align_context.moveTo(x + i, y);
        align_context.lineTo(x + i, y + (height / 2));
        align_context.stroke();
    }

}


//
//
//
function visualizeMultiple(order) {
    // Visual Constants
    var divMargin = 10; // Left margin in div block in pixels
    var alignmentGlyphHeight = 1;
    var alignmentSpacing = alignmentGlyphHeight + 1;
    var rulerHeight = 8;
    var rulerVerticalMargin = 10;
    
    var viewWidth = align_canvas.width - divMargin; // Width of reference sequence in pixels
    var xScale = viewWidth / data.length;
    var alignments = data.alignments;
    var qualWidthBP = data.qualityBlockLen;
    
    // Clear overlayed canvases
    align_context.clearRect(0, 0, align_canvas.width,
    align_canvas.height);
    guide_context.clearRect(0, 0, guide_canvas.width,
    guide_canvas.height);

    // Select ordering
    if (order == "orient") {
        alignments.sort(function (a, b) {
            if (a[4] === b[4]) {
                if (a[4] === "R") {
                    if (a[1] === b[1]) {
                        return (b[2] - a[2]);
                    } else {
                        return (a[1] - b[1]);
                    }
                } else {
                    if (a[1] === b[1]) {
                        return (a[2] - b[2]);
                    } else {
                        return (b[1] - a[1]);
                    }
                }
            } else {
                return a[4] < b[4] ? -1 : a[4] > b[4] ? 1 : 0;
            }
        });
    } else if ( order == "end" )
    {
        alignments.sort(function (a, b) {
            if ((a[1] + a[2]) == (b[1] + b[2])) {
                return (a[1]  - b[1]);
            } else {
                return ((b[1] + b[2]) - (a[1] + a[2]));
            }
        });        
    } else if ( order == "div" )
    {
        alignments.sort(function (a, b) {
            return (a[5] - b[5]);
        });        
     }else {
        alignments.sort(function (a, b) {
            if (a[1] === b[1]) {
                return (b[2] - a[2]);
            } else {
                return (a[1] - b[1]);
            }
        });
    }


    var curY = 0;
    var referenceDrawn = 0;
    if (order != "orient") {
        ruler(divMargin, 0, viewWidth, rulerHeight, 1, 946, 10, 100);
        currRulerY = 0;
        curY = rulerVerticalMargin + rulerHeight;        
        referenceDrawn = 1;
    }
    for (var i = 0; i < alignments.length; i += 1) {
        if (referenceDrawn == 0 && alignments[i][4] == "R") {
            curY = curY + rulerVerticalMargin;
            ruler(divMargin, curY + (i * alignmentSpacing),
            viewWidth, rulerHeight, 1, 946, 10, 100);
            currRulerY = curY + (i * alignmentSpacing);
            curY = curY + rulerHeight + rulerVerticalMargin;
            referenceDrawn = 1;
        }

        var xOffset = alignments[i][1];
        var qualities = alignments[i][3];
        var qualIdx = 0;
        for (var j = 0; j < alignments[i][2]; j += qualWidthBP) {

            // TODO fix this indexing error
            if (qualIdx < qualities.length) {

                var grd = align_context.createLinearGradient(
                            0, 0, (xScale * qualWidthBP), 0 );

                var grad = qualColor[qualities[qualIdx] - 1] + " - ";
                grd.addColorStop(0, qualColor[qualities[qualIdx] - 1]);
                if (qualIdx == qualities.length - 1) {
                    grd.addColorStop(1, 
                                   qualColor[qualities[qualIdx] - 1]);
                    grad = grad + qualColor[qualities[qualIdx] - 1];
                } else {
                    grd.addColorStop(1, 
                                 qualColor[qualities[qualIdx + 1] - 1]);
                    grad = grad + qualColor[qualities[qualIdx + 1] - 1];
                }
                align_context.fillStyle = grd;
                align_context.fillRect(divMargin + (xOffset * xScale) +
                                       (j * xScale),
                curY + (i * alignmentSpacing), (xScale * qualWidthBP),
                alignmentGlyphHeight);
            }

            qualIdx++;
        }
    }
}


//
//  Show alignment
//
visualizeMultiple("norm");
';
  return ( $javascript );
}

sub inlineTableJS
{
  my $javascript = ' 
// ===================================================================
// Author: Denis Howlett <feedback@isocra.com>
// WWW: http://www.isocra.com/
// ===================================================================
                
/** Keep hold of the current table being dragged */
var currenttable = null;
            
/** Capture the onmousemove so that we can see if a row from the current
 *  table if any is being dragged.
 * @param ev the event (for Firefox and Safari, otherwise we use window.event for IE)
 */ 
document.onmousemove = function(ev){
    if (currenttable && currenttable.dragObject) {
        ev   = ev || window.event;
        var mousePos = currenttable.mouseCoords(ev);
        var y = mousePos.y - currenttable.mouseOffset.y;
        if (y != currenttable.oldY) {
            // work out if were going up or down...
            var movingDown = y > currenttable.oldY;
            // update the old value
            currenttable.oldY = y;
            // update the style to show were dragging
            currenttable.dragObject.style.backgroundColor = "#eee";
            // If were over a row then move the dragged row to there so that the user sees the
            // effect dynamically
            var currentRow = currenttable.findDropTargetRow(y);
            if (currentRow) {
                if (movingDown && currenttable.dragObject != currentRow) {
                    currenttable.dragObject.parentNode.insertBefore(currenttable.dragObject, currentRow.nextSibling);
                } else if (! movingDown && currenttable.dragObject != currentRow) {
                    currenttable.dragObject.parentNode.insertBefore(currenttable.dragObject, currentRow);
                }
            }
        }

        return false;
    }
}

// Similarly for the mouseup
document.onmouseup   = function(ev){
    if (currenttable && currenttable.dragObject) {
        var droppedRow = currenttable.dragObject;
        // If we have a dragObject, then we need to release it,
        // The row will already have been moved to the right place so we just reset stuff
        droppedRow.style.backgroundColor = "transparent";
        currenttable.dragObject   = null;
        // And then call the onDrop method in case anyone wants to do any post processing
        currenttable.onDrop(currenttable.table, droppedRow);
        currenttable = null; // let go of the table too
    }
}


/** get the source element from an event in a way that works for IE and Firefox and Safari
 * @param evt the source event for Firefox (but not IE--IE uses window.event) */
function getEventSource(evt) {
    if (window.event) {
        evt = window.event; // For IE
        return evt.srcElement;
    } else {
        return evt.target; // For Firefox
    }
}

/**
 * Encapsulate table Drag and Drop in a class. Well have this as a Singleton
 * so we dont get scoping problems.
 */
function TableDnD() 
{
    /** Keep hold of the current drag object if any */
    this.dragObject = null;
    /** The current mouse offset */
    this.mouseOffset = null;
    /** The current table */
    this.table = null;
    /** Remember the old value of Y so that we dont do too much processing */
    this.oldY = 0;

    /** Initialise the drag and drop by capturing mouse move events */
    this.init = function(table) {
        this.table = table;
        var rows = table.tBodies[0].rows; //getElementsByTagName("tr")
        for (var i=0; i<rows.length; i++) {
                        // John Tarr: added to ignore rows that Ive added the NoDnD attribute to (Category and Header rows)
                        var nodrag = rows[i].getAttribute("NoDrag")
                        if (nodrag == null || nodrag == "undefined") { //There is no NoDnD attribute on rows I want to drag
                                this.makeDraggable(rows[i]);
                        }
        }
    }

    /** This function is called when you drop a row, so redefine it in your code
        to do whatever you want, for example use Ajax to update the server */
    this.onDrop = function(table, droppedRow) {
        // Do nothing for now
    }

        /** Get the position of an element by going up the DOM tree and adding up all the offsets */
    this.getPosition = function(e){
        var left = 0;
        var top  = 0;
                /** Safari fix -- thanks to Luis Chato for this! */
                if (e.offsetHeight == 0) {
                        /** Safari 2 doesnt correctly grab the offsetTop of a table row
                            this is detailed here:
                            http://jacob.peargrove.com/blog/2006/technical/table-row-offsettop-bug-in-safari/
                            the solution is likewise noted there, grab the offset of a table cell in the row - the firstChild.
                            note that firefox will return a text node as a first child, so designing a more thorough
                            solution may need to take that into account, for now this seems to work in firefox, safari, ie */
                        e = e.firstChild; // a table cell
                }
        while (e.offsetParent){
            left += e.offsetLeft;
            top  += e.offsetTop;
            e     = e.offsetParent;
        }

        left += e.offsetLeft;
        top  += e.offsetTop;

        return {x:left, y:top};
    }

        /** Get the mouse coordinates from the event (allowing for browser differences) */
    this.mouseCoords = function(ev){
        if(ev.pageX || ev.pageY){
            return {x:ev.pageX, y:ev.pageY};
        }
        return {
            x:ev.clientX + document.body.scrollLeft - document.body.clientLeft,
            y:ev.clientY + document.body.scrollTop  - document.body.clientTop
        };
    }


        /** Given a target element and a mouse event, get the mouse offset from that element.
                To do this we need the elements position and the mouse position */
    this.getMouseOffset = function(target, ev){
        ev = ev || window.event;

        var docPos    = this.getPosition(target);
        var mousePos  = this.mouseCoords(ev);
        return {x:mousePos.x - docPos.x, y:mousePos.y - docPos.y};
    }



        /** Take an item and add an onmousedown method so that we can make it draggable */
    this.makeDraggable = function(item) {
        if(!item) return;
        var self = this; // Keep the context of the TableDnd inside the function
        item.onmousedown = function(ev) {
            // Need to check to see if we are an input or not, if we are an input, then
            // return true to allow normal processing
            var target = getEventSource(ev);
            if (target.tagName == "INPUT" || target.tagName == "SELECT") return true;
            currenttable = self;
            self.dragObject  = this;
            self.mouseOffset = self.getMouseOffset(this, ev);
            return false;
        }
        item.style.cursor = "move";
    }

    /** Were only worried about the y position really, because we can only move rows up and down */
    this.findDropTargetRow = function(y) {
        var rows = this.table.tBodies[0].rows;
                for (var i=0; i<rows.length; i++) {
                        var row = rows[i];
                        // John Tarr added to ignore rows that Ive added the NoDnD attribute to (Header rows)
                        var nodrop = row.getAttribute("NoDrop");
                        if (nodrop == null || nodrop == "undefined") {  //There is no NoDnD attribute on rows I want to drag
                                var rowY    = this.getPosition(row).y;
                                var rowHeight = parseInt(row.offsetHeight)/2;
                                if (row.offsetHeight == 0) {
                                        rowY = this.getPosition(row.firstChild).y;
                                        rowHeight = parseInt(row.firstChild.offsetHeight)/2;
                                }
                                // Because we always have to insert before, we need to offset the height a bit
                                if ((y > rowY - rowHeight) && (y < (rowY + rowHeight))) {
                                        // thats the row were over
                                        return row;
                                }
                        }
                }
                return null;
        }
}';
  return $javascript;
}

sub getMSASummaryViewJSON
{
  my $object = shift;

  my $maxNameLen      = 0;
  my $qualityBlockLen = 10;

  # Get ungapped reference sequence length
  my $refSeq      = $object->getReferenceSeq();
  my $refSeqNoIns = $refSeq;
  $refSeqNoIns =~ s/-//g;
  my $refLen = length( $refSeqNoIns );

  # Build datastructure so we can export it in javascript syntax
  my %data = ();

  $object->kimuraDivergence();
  $data{'num_alignments'}  = $object->getNumAlignedSeqs();
  $data{'length'}          = $refLen;
  $data{'qualityBlockLen'} = $qualityBlockLen;
  $data{'alignments'}      = [];

  for ( my $i = 0 ; $i < $object->getNumAlignedSeqs() ; $i++ )
  {

    # Count referene up to aligned start
    my $alignedStart = $object->getAlignedStart( $i );
    my $refStart     = 0;
    for ( my $j = 0 ; $j <= $alignedStart ; $j++ )
    {
      if ( substr( $refSeq, $j, 1 ) ne "-" )
      {
        $refStart++;
      }
    }
    my $seq     = $object->getAlignedSeq( $i );
    my $ins     = 0;
    my $del     = 0;
    my $mut     = 0;
    my $aligned = 0;
    my @scores  = ();
    for ( my $j = 0 ; $j < length( $seq ) ; $j++ )
    {
      my $rChar = substr( $refSeq, $j + $alignedStart, 1 );
      my $aChar = substr( $seq, $j, 1 );
      if ( $rChar eq "-" )
      {
        $ins++ if ( $aChar ne "-" );
        next;
      }
      $aligned++;

      # TODO: Consider the order of these two conditionals.  It seems that del
      # would never be invoked.
      if ( $rChar ne $aChar )
      {
        $mut++;
      } elsif ( $aChar eq "-" )
      {
        $del++;
      }
      if ( $aligned % $qualityBlockLen == 0 )
      {
        my $score = $qualityBlockLen - $mut + $del;
        $score-- if ( $ins );
        $score = 1 if ( $score < 1 );
        push @scores, $score;
        $ins = 0;
        $del = 0;
        $mut = 0;
      }
    }
    if ( $mut || $del || $ins )
    {
      my $score = $qualityBlockLen - $mut + $del;
      $score-- if ( $ins );
      $score = 1 if ( $score < 1 );
      push @scores, $score;
    }

    my $name   = $object->getAlignedName( $i );
    my $div    = sprintf( "%0.2f", $object->getAlignedDiv( $i ) );
    my $start  = $object->getAlignedSeqStart( $i );
    my $end    = $object->getAlignedSeqEnd( $i );
    my $orient = "F";
    if ( $object->getAlignedOrientation( $i ) eq "-" )
    {
      $orient = "R";
    }
    push @{ $data{'alignments'} },
        [ $name, $refStart, $aligned, [ @scores ], $orient, $div, $start,
          $end ];

  }

  my $json = JSON::PP->new->utf8;
  return ( $json->encode( \%data ) );
}

sub generateJavascriptSummaryAndAlignmentViewer
{
  my $mAlign = shift;
  my $MOUT   = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $OUT = *STDOUT;
  if ( defined $MOUT )
  {
    if ( ref( $MOUT ) !~ /GLOB|FileHandle/ )
    {
      print "::$subroutine(" . ") Opening file " . $MOUT . "\n"
          if ( $DEBUG );
      open $OUT, $MOUT
          or die
          . "$subroutine("
          . ": Unable to open "
          . "results file: $MOUT : $!";
    } else
    {
      $OUT = $MOUT;
    }
  }

  my $json = JSON::PP->new->utf8;

  # Start building detail json data
  my %detailData = ();
  $detailData{'alignment'}  = ();
  $detailData{'alignWidth'} = length( $mAlign->getReferenceSeq() );

  # Save the alignment reference
  push @{ $detailData{'alignment'} },
      {
        'id'       => 'reference',
        'sequence' => $mAlign->getReferenceSeq()
      };

  # Save the score data
  my ( $columns, $scoreArray ) = $mAlign->getLowScoringAlignmentColumns();
  $detailData{'alignmentScore'} = ();
  for ( my $j = 0 ; $j <= $#{$scoreArray} ; $j++ )
  {
    my $num = int( $scoreArray->[ $j ] * 10 ) / 10;
    push @{ $detailData{'alignmentScore'} }, $num;
    $scoreArray->[ $j ] = $num;
  }

  # Save the actual alignments
  for ( my $i = 0 ; $i < $mAlign->getNumAlignedSeqs ; $i++ )
  {
    my $seq   = $mAlign->getAlignedSeq( $i );
    my $name  = $mAlign->getAlignedName( $i );
    my $start = $mAlign->getAlignedStart( $i );
    push @{ $detailData{'alignment'} },
        {
          'id'       => $name,
          'start'    => $start,
          'sequence' => $seq
        };
  }

  # Start building summary json data
  my %summaryData = ();

  # Calculate and set divergence values in MultAlign
  $mAlign->kimuraDivergence();

  # Get ungapped reference sequence length
  my $refSeq      = $mAlign->getReferenceSeq();
  my $refSeqNoIns = $refSeq;
  $refSeqNoIns =~ s/-//g;
  my $refLen          = length( $refSeqNoIns );
  my $qualityBlockLen = 10;
  $summaryData{'num_alignments'}  = $mAlign->getNumAlignedSeqs();
  $summaryData{'length'}          = $refLen;
  $summaryData{'qualityBlockLen'} = $qualityBlockLen;
  $summaryData{'alignments'}      = [];
  for ( my $i = 0 ; $i < $mAlign->getNumAlignedSeqs ; $i++ )
  {

    # Count referene up to aligned start
    my $alignedStart = $mAlign->getAlignedStart( $i );
    my $refStart     = 0;
    for ( my $j = 0 ; $j <= $alignedStart ; $j++ )
    {
      if ( substr( $refSeq, $j, 1 ) ne "-" )
      {
        $refStart++;
      }
    }
    my $seq     = $mAlign->getAlignedSeq( $i );
    my $ins     = 0;
    my $del     = 0;
    my $mut     = 0;
    my $aligned = 0;
    my @scores  = ();
    for ( my $j = 0 ; $j < length( $seq ) ; $j++ )
    {
      my $rChar = substr( $refSeq, $j + $alignedStart, 1 );
      my $aChar = substr( $seq, $j, 1 );
      if ( $rChar eq "-" )
      {
        $ins++ if ( $aChar ne "-" );
        next;
      }
      $aligned++;
      if ( $rChar ne $aChar )
      {
        $mut++;
      } elsif ( $aChar eq "-" )
      {
        $del++;
      }
      if ( $aligned % $qualityBlockLen == 0 )
      {
        my $score = $qualityBlockLen - $mut + $del;
        $score-- if ( $ins );
        $score = 1 if ( $score < 1 );
        push @scores, $score;
        $ins = 0;
        $del = 0;
        $mut = 0;
      }
    }
    if ( $mut || $del || $ins )
    {
      my $score = $qualityBlockLen - $mut + $del;
      $score-- if ( $ins );
      $score = 1 if ( $score < 1 );
      push @scores, $score;
    }

    my $name   = $mAlign->getAlignedName( $i );
    my $div    = sprintf( "%0.2f", $mAlign->getAlignedDiv( $i ) );
    my $start  = $mAlign->getAlignedSeqStart( $i );
    my $end    = $mAlign->getAlignedSeqEnd( $i );
    my $orient = "F";
    if ( $mAlign->getAlignedOrientation( $i ) eq "-" )
    {
      $orient = "R";
    }
    push @{ $summaryData{'alignments'} },
        [ $name, $refStart, $aligned, [ @scores ], $orient, $div, $start,
          $end ];

  }

  # Begin writing the HTML
  print $OUT "<html>
<H1>Summary View</H1>
<button onClick=\"mySummary.render('orient');\">Orientation Sort</button>
<button onClick=\"mySummary.render('norm');\">Normal Sort</button>
<button onClick=\"mySummary.render('end');\">End Sort</button>
<button onClick=\"mySummary.render('div');\">Divergence Sort</button>
<p>    
<div id=\"canvasesdiv\" style=\"position:relative\">
  <canvas id=\"alignment_canvas\" width=\"800\" height=\"1600\" style=\"z-index:1;position:absolute;left:0px;top:0px;\">Canvas not supported</canvas>
  <canvas id=\"detail_canvas\" width=\"800\" height=\"1600\" style=\"z-index:2;position:absolute;left:0px;top:0px;\">Canvas not supported</canvas>
  <canvas id=\"guideline_canvas\" width=\"800\" height=\"1600\" style=\"z-index:3;position:absolute;left:0px;top:0px;\">Canvas not supported</canvas>
</div>
<p>
<h1>Detail View</h1>
<button onClick=\"myViewer.setViewType('norm');\">Normal View</button>
<button onClick=\"myViewer.setViewType('diffs');\">Difference View</button>
<canvas id=\"canvas\" width=\"1500\" height=\"700\"></canvas>
<p><script>\n";

  # Inline javascript dependencies
  open IN, "<$FindBin::RealBin/javascript/zynga-1.2.2-10/Animate.js"
      or die "Could not inline $FindBin::RealBin/javascript/"
      . "zynga-1.2.2-10/Animate.js file!";
  while ( <IN> )
  {
    print $OUT "$_";
  }
  close IN;
  open IN, "<$FindBin::RealBin/javascript/zynga-1.2.2-10/Scroller.js"
      or die "Could not inline $FindBin::RealBin/javascript/"
      . "zynga-1.2.2-10/Scroller.js file!";
  while ( <IN> )
  {
    print $OUT "$_";
  }
  close IN;

  print $OUT "var summaryData = ";
  my $jsonStr = $json->encode( \%summaryData );

  # For using jsfiddle
  #$jsonStr =~ s/\},/\},\n/g;
  #$jsonStr =~ s/\[/\n[/g;
  print $OUT "$jsonStr;\n";

  open IN, "<$FindBin::RealBin/javascript/isb/AlignmentSummary.js"
      or die "Could not inline $FindBin::RealBin/javascript/"
      . "isb/AlignmentSummary.js file!";
  while ( <IN> )
  {
    print $OUT "$_";
  }
  close IN;

  print $OUT "\n\n";

  print $OUT "var detailData = ";
  $jsonStr = $json->encode( \%detailData );

  # For using jsfiddle
  #$jsonStr =~ s/\},/\},\n/g;
  print $OUT "$jsonStr;\n";

  open IN, "<$FindBin::RealBin/javascript/isb/AlignmentViewer.js"
      or die "Could not inline $FindBin::RealBin/javascript/"
      . "isb/AlignmentViewer.js file!";
  while ( <IN> )
  {
    print $OUT "$_";
  }
  close IN;

  print $OUT "var mySummary = new AlignmentSummary( "
      . "  document.getElementById('alignment_canvas'), "
      . "  document.getElementById('guideline_canvas'), "
      . "  document.getElementById('detail_canvas'), "
      . " summaryData, {});\n";

  print $OUT "var myViewer = new AlignmentViewer( "
      . "document.getElementById('canvas'), detailData, {} );\n";
  print $OUT "</script></html>";
}

1;
