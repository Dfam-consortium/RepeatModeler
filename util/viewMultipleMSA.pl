#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) viewMultipleMSA
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
###############################################################################

=head1 NAME

viewMultipleMSA - View multiple MSA alignments using HTML 

=head1 SYNOPSIS

  viewMultipleMSA.pl msa.fa [msa.fa ..]

=head1 DESCRIPTION

  No description yet....

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2020 Robert Hubley, Institute for Systems Biology

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
use MultAln;

# RepeatMasker Libraries
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
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
my $Version = $RepModelConfig::VERSION;

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
                    '-fullmsa',
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

my $MOUT;
open $MOUT,">MultMSA.html" or die;

my $idx = 1;
foreach my $file ( @ARGV ) {
  print "Processing $file\n";
  my $malign = readMSA($file);
  generateJavascriptSummaryAndAlignmentViewer( $malign, $MOUT, $file, $idx );
  $idx++;
}

sub readMSA {
  my $file = shift;
  my $seq;
  my $id;
  my @seqs;
  open my $IN, "<$file" or die "Could not open $file for reading";
  # Simple FASTA reader
  my %idHash = ();
  while (<$IN>) {
    if ( /^>(\S+)/ )
    {
      my $tmpID = $1;
      if ( defined $idHash{$tmpID} ) {
        my $ver = 1;
        while ( defined $idHash{$tmpID . "_$ver"} )
        {
          $ver++;
        }
        warn "WARN File contains a duplicate identifier \"$tmpID\".  A suffix of \"_$ver\"\n" .
             "     will be appended to this occurence.\n";
        $tmpID = $tmpID . "_$ver";
      }
      $idHash{$tmpID}++;
      if ( $seq )
      {
        $seq = uc($seq);
        # Convert prefix/suffix "-"s to spacers
        if ( $seq =~ /^(\-+)/ ){
          substr($seq,0,length($1)) = " "x(length($1));
        }
        if ( $seq =~ /(\-+)$/ ) {
          substr($seq,length($seq)-length($1)-1) = " "x(length($1));
        }
#print "$id : >$seq<\n";
        push @seqs, [ $id, $seq ];
      }
      $seq = "";
      $id = $tmpID;
      next;
    }
    s/[\s\n\r]+//g;
    $seq .= $_;
  }
  if ( $seq )
  {
    # Convert prefix/suffix "-"s to spacers
    if ( $seq =~ /^(\-+)/ ){
      substr($seq,0,length($1)) = " "x(length($1));
    }
    if ( $seq =~ /(\-+)$/ ) {
      substr($seq,length($seq)-length($1)-1) = " "x(length($1));
    }

    $seq = uc($seq);
    push @seqs, [ $id, $seq ];
  }
  close $IN;
  my $malign = MultAln->new( sequences => \@seqs );
  return $malign;
}




######################## S U B R O U T I N E S ############################

sub generateJavascriptSummaryAndAlignmentViewer {
  my $mAlign = shift;
  my $MOUT   = shift;
  my $filename = shift;
  my $idx = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $OUT = *STDOUT;
  if ( defined $MOUT ) {
    if ( ref( $MOUT ) !~ /GLOB|FileHandle/ ) {
      print "::$subroutine(" . ") Opening file " . $MOUT . "\n"
          if ( $DEBUG );
      open $OUT, $MOUT
          or die
          . "$subroutine("
          . ": Unable to open "
          . "results file: $MOUT : $!";
    }
    else {
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
  for ( my $j = 0 ; $j <= $#{$scoreArray} ; $j++ ) {
    my $num = int( $scoreArray->[ $j ] * 10 ) / 10;
    push @{ $detailData{'alignmentScore'} }, $num;
    $scoreArray->[ $j ] = $num;
  }

  # Save the actual alignments
  for ( my $i = 0 ; $i < $mAlign->getNumAlignedSeqs ; $i++ ) {
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
  if ( $options{'fullmsa'} ) {
    $refLen = length($refSeq);
  }
  my $qualityBlockLen = 10;
  $summaryData{'num_alignments'}  = $mAlign->getNumAlignedSeqs();
  $summaryData{'length'}          = $refLen;
  $summaryData{'qualityBlockLen'} = $qualityBlockLen;
  $summaryData{'alignments'}      = [];
  for ( my $i = 0 ; $i < $mAlign->getNumAlignedSeqs ; $i++ ) {

    # Count referene up to aligned start
    my $alignedStart = $mAlign->getAlignedStart( $i );
    my $refStart     = 0;
    if ( $options{'fullmsa'} ) {
      $refStart = $alignedStart;
    }else{
      for ( my $j = 0 ; $j <= $alignedStart ; $j++ ) {
        if ( substr( $refSeq, $j, 1 ) ne "-" ) {
          $refStart++;
        }
      }
    }
    my $seq     = $mAlign->getAlignedSeq( $i );
    my $ins     = 0; # Reference Insertion Count: e.g Ref=-, Seq=A
    my $del     = 0; # Reference Deletion Count: e.g Ref=A, Seq=-
    my $mut     = 0; # Reference Mutation Count: e.g Ref=A, Seq=C
    my $una     = 0; # Reference Unaligned Count: e.g Ref=-, Seq=-
    my $idt     = 0; # Reference Identity Count: e.g Ref=A, Seq=A
    my $currBlockLen = 0;
    my $totalBlockLen = 0;
    my @scores  = ();
    for ( my $j = 0 ; $j < length( $seq ) ; $j++ ) {
      my $rChar = substr( $refSeq, $j + $alignedStart, 1 );
      my $aChar = substr( $seq, $j, 1 );
      
      if ( $rChar eq "-" ) {
        $currBlockLen++ if ( $options{'fullmsa'} );
        if ( $aChar eq "-" ) {
          $una++;
        }else {
          $ins++;
        }
      }else {
        $currBlockLen++;
        if ( $aChar eq "-" ) {
          $del++;
        }elsif ( $aChar ne $rChar ) {
          $mut++;
        }else {
          $idt++;
        }
      }
      if ( $currBlockLen == $qualityBlockLen ) {
        $totalBlockLen += $currBlockLen;
        my $score;
        if ( $options{'fullmsa'} && $una == $qualityBlockLen ) {
          # Case where we want to see insertions *and* ref=-, seq=- for the entire
          # block.  In this case it should just be transparent.
          $score = 0;
        }else {
          #$score = $qualityBlockLen - ($mut + $del);
          $score = $idt;
          $score-- if ( $ins );
          $score = 1 if ( $score < 1 );
        }
        push @scores, $score;
        $ins = 0;
        $del = 0;
        $mut = 0;
        $una = 0;
        $idt = 0;
        $currBlockLen = 0;
      }
    }
    if ( $mut || $del || $ins || $idt || $una ) {
      $totalBlockLen += $currBlockLen;
      my $score;
      if ( $options{'fullmsa'} && $una == $qualityBlockLen ) {
        # Case where we want to see insertions *and* ref=-, seq=- for the entire
        # block.  In this case it should just be transparent.
        $score = 0;
      }else {
        #$score = $qualityBlockLen - ($mut + $del);
        $score = $idt;
        $score-- if ( $ins );
        $score = 1 if ( $score < 1 );
      }
      push @scores, $score;
    }

    my $name   = $mAlign->getAlignedName( $i );
    my $div    = sprintf( "%0.2f", $mAlign->getAlignedDiv( $i ) );
    my $start  = $mAlign->getAlignedSeqStart( $i );
    my $end    = $mAlign->getAlignedSeqEnd( $i );
    my $orient = "F";
    if ( $mAlign->getAlignedOrientation( $i ) eq "-" ) {
      $orient = "R";
    }
    push @{ $summaryData{'alignments'} },
        [ $name, $refStart, $totalBlockLen, [ @scores ], $orient, $div, $start,
          $end ];

  }

  # Begin writing the HTML
  print $OUT "<html>\n";
  if ( $options{'fullmsa'} ) {
    print $OUT "<H1>Full MSA View: $filename</H1>\n";
  }else {
    print $OUT "<H1>Consensus MSA View: $filename</H1>\n";
  }
print $OUT "
<button onClick=\"mySummary_$idx.render('orient');\">Orientation Sort</button>
<button onClick=\"mySummary_$idx.render('norm');\">Normal Sort</button>
<button onClick=\"mySummary_$idx.render('end');\">End Sort</button>
<button onClick=\"mySummary_$idx.render('div');\">Divergence Sort</button>
<p>    
<div id=\"canvasesdiv\" style=\"position:relative\">
  <canvas id=\"alignment_canvas_$idx\" width=\"800\" height=\"1600\" style=\"z-index:1;position:absolute;left:0px;top:0px;\">Canvas not supported</canvas>
  <canvas id=\"detail_canvas_$idx\" width=\"800\" height=\"1600\" style=\"z-index:2;position:absolute;left:0px;top:0px;\">Canvas not supported</canvas>
  <canvas id=\"guideline_canvas_$idx\" width=\"800\" height=\"1600\" style=\"z-index:3;position:absolute;left:0px;top:0px;\">Canvas not supported</canvas>
</div>
<p>\n";


print $OUT "<p><script>\n";

  print $OUT "var summaryData_$idx = ";
  my $jsonStr = $json->encode( \%summaryData );

  # For using jsfiddle
  #$jsonStr =~ s/\},/\},\n/g;
  #$jsonStr =~ s/\[/\n[/g;
  print $OUT "$jsonStr;\n";

  open IN, "<$FindBin::RealBin/javascript/isb/AlignmentSummary.js"
      or die "Could not inline $FindBin::RealBin/javascript/"
      . "isb/AlignmentSummary.js file!";
  while ( <IN> ) {
    print $OUT "$_";
  }
  close IN;

  print $OUT "\n\n";

  print $OUT "var mySummary_$idx = new AlignmentSummary( "
      . "  document.getElementById('alignment_canvas_$idx'), "
      . "  document.getElementById('guideline_canvas_$idx'), "
      . "  document.getElementById('detail_canvas_$idx'), "
      . " summaryData_$idx, {});\n";

#  print $OUT "var myViewer = new AlignmentViewer( "
#      . "document.getElementById('canvas'), detailData, {} );\n";
  print $OUT "</script></html>";
}

1;
