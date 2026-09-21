#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) RepeatUtil.pm
##  Authors:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Arian Smit         asmit@systemsbiology.org
##  Description:
##      A module containing several useful subroutines
##      used by the RepeatModeler suite of programs.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2004 Developed by
#* Robert Hubley, Arian Smit and Arnie Kas.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
#  ChangeLog:
#
#    $Log: RepeatUtil.pm,v $
#    Revision 1.29  2017/04/05 00:03:31  rhubley
#    Cleanup before a distribution
#
#
###############################################################################

=head1 NAME

RepeatUtil.pm - Library functions for RepeatModeler

=head1 SYNOPSIS

use RepeatUtil;

Usage:


=head1 DESCRIPTION

=head1 SEE ALSO

=over 4

RepeatModeler 

=back

=head1 COPYRIGHT

Copyright 2005 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>
Arian Smit <asmit@systemsbiology.org>

=head1 ATTRIBUTES

=cut

#
# Module Dependence
#
package RepeatUtil;
use strict;
use FindBin;
use lib $FindBin::RealBin;
use Data::Dumper;
use Carp;
use File::Basename;

# RepeatMasker Libraries
use RepModelConfig;
use CrossmatchSearchEngine;
use NCBIBlastSearchEngine;
use MultAln;
use SeedAlignment;
use SeqDBI;
use SearchResultCollection;
use SimpleBatcher;
use ThreadedTaskSimple;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "RepeatUtil";
my $DEBUG = 0;
$DEBUG = 1 if ( $RepModelConfig::DEBUGALL == 1 );
my $config        = $RepModelConfig::configuration;


##---------------------------------------------------------------------##

=head2 ncbiMaskDatabaseNativeMT()

  Use: ncbiMaskDatabaseNativeMT( 
                     rmblastnPath => "/usr/local/rmblast/rmblastn",
                     fastaFile => "/jo/bob/seq.fa",
                     consensi => "/jo/bob/lib/reps.fa",
                     workingDir => "/jo/bob/round-3",
                     threads => 4,
                     [instSeqFile => "/jo/bob/instances.out"],
                      );

  Mask the fastaFile (as the query) using the TE consensi (as the
  database) with rmblastn's query threading (mt_mode = 1), which
  RMBlast has had since 2.13.  prepareSubject() prepares the consensi in
  workingDir, whatever that means for the installed rmblastn series.

  This routine instantiates its own searchEngine in order to manage
  its settings more effectively.

=cut

##---------------------------------------------------------------------##
sub ncbiMaskDatabaseNativeMT {
  my %parameters = @_;

  my $fName = "ncbiMaskDatabaseNativeMT";

  # Parameter checking
  die $CLASS
      . "::$fName(): Missing or invalid rmblastnPath "
      . "parameter!\n"
      if (    !defined $parameters{'rmblastnPath'}
           || !-x $parameters{'rmblastnPath'} );
  my $rmblastnPath = $parameters{'rmblastnPath'};

  die $CLASS . "::$fName(): Missing workingDir parameter!\n"
      if (    !defined $parameters{'workingDir'}
           || !-d $parameters{'workingDir'} );
  my $workingDir = $parameters{'workingDir'};

  die $CLASS . "::$fName(): Missing fastaFile parameter!\n"
      if (    !defined $parameters{'fastaFile'}
           || !-s $parameters{'fastaFile'} );
  my $fastaFile = $parameters{'fastaFile'};

  die $CLASS . "::$fName(): Missing consensi parameter!\n"
      if (    !defined $parameters{'consensi'}
           || !-s $parameters{'consensi'} );
  my $consensi = $parameters{'consensi'};

  die $CLASS . "::$fName(): Missing threads parameter!\n"
      if (    !defined $parameters{'threads'} );
  my $threads = $parameters{'threads'};

  my $searchEngine =
        NCBIBlastSearchEngine->new( pathToEngine => $rmblastnPath );

  if ( $searchEngine->getVersion() =~ /.*(\d+)\.(\d+)\.\d+/ ){
    unless ( $1 > 2 || ($1 == 2 && $2 >= 13) )
    {
      die $CLASS . "::$fName(): rmblast version is incompatible with this function!\n";
    }
  }else {
    die $CLASS . "::$fName(): could not obtain rmblast version!\n";
  }

  ## TODO FORMALIZE once evaluated
  ## 2025-05-19 : Lowering to improve masking performance
  #$searchEngine->setMinScore( 250 );
  $searchEngine->setMinScore( 200 );
  # TODO one for each
  $searchEngine->setMinScore($parameters{'minScore'}) 
     if ( exists $parameters{'minScore'} );
  
  $searchEngine->setGenerateAlignments( 0 );
  $searchEngine->setGapInit( -25 );
  $searchEngine->setBandwidth( 10 );    # Changes gapW=31
  $searchEngine->setInsGapExt( -5 );
  $searchEngine->setDelGapExt( -5 );
  $searchEngine->setMinMatch( 7 );
  $searchEngine->setScoreMode( SearchEngineI::complexityAdjustedScoreMode );
  $searchEngine->setMaskLevel( 90 );
  $searchEngine->setCores($threads);
  $searchEngine->setThreadByQuery(1);
  $searchEngine->setTempDir($workingDir);
  $searchEngine->setMatrix(
                     "$FindBin::RealBin/Matrices/ncbi/nt/comparison.matrix" );

  # Setup the temporary database.  The consensi file grows between
  # rounds, so force a rebuild rather than reuse an index left by an
  # earlier call.
  my $subject = $searchEngine->prepareSubject( $consensi,
                                               outputDir   => $workingDir,
                                               dbName      => "tmpConsDB",
                                               parseSeqIDs => 1,
                                               dbVersion   => 4,
                                               force       => 1 );

  my $maskDB = FastaDB->new( fileName => $fastaFile,
                             openMode => SeqDBI::ReadOnly );
  my %maskSeqs   = ();
  my @maskSeqIDs = $maskDB->getIDs();
  foreach my $seqID ( @maskSeqIDs ) {
    my $seq  = $maskDB->getSequence( $seqID );
    $maskSeqs{$seqID} = $seq; 
  }
  undef $maskDB;
 
  # Setup the Query/Subject
  $searchEngine->setQuery($fastaFile);
  $searchEngine->setSubject($subject);

  my $INST;
  if ( defined $parameters{'instSeqFile'} ) {
    open $INST, ">$parameters{'instSeqFile'}";
  }

  my $repeatsMasked = 0;
  my $totalMasked = 0;

  ## DEBUG
  ##print "Running: " . $searchEngine->getParameters() . "\n";
  ## DEBUG
  my ( $status, $resultCollection ) = $searchEngine->search();
  if ( $status ) {
    print STDERR "\nERROR from search engine (", $? >> 8, ") \n";
  }
  elsif ( $resultCollection->size() > 0 ) {
    print "    -- Collecting " . $resultCollection->size() . " ranges...\n";
    $resultCollection->sort(
            sub ($$) {
               ($_[ 0 ]->getQueryName() cmp $_[ 1 ]->getQueryName()) ||
               ($_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart()) ||
               ($_[ 1 ]->getQueryEnd() <=> $_[ 0 ]->getQueryEnd());
                     });
    my $prevStart = -1;
    my $prevEnd = -1;
    my $prevID = "";
    for ( my $k = 0 ; $k < $resultCollection->size(); $k++ ) {
      my $res = $resultCollection->get( $k ); 
      my $seqID = $res->getQueryName();
      my $start = $res->getQueryStart();
      my $end = $res->getQueryEnd();
      my $len = $end - $start + 1;

      # Reset previous stats if we cross a sequence boundary
      if ( $seqID ne $prevID ) {
        $prevID = "";
        $prevStart = -1;
        $prevEnd = -1;
      }

      # Cacluate actual sequence masked (accounting for overlaps)
      my $actual_masked = $len;
      if ( $prevStart > 0 ) {
        if ( $prevEnd > $start ) {
          if ( $prevEnd > $end ) {
            $actual_masked = 0;
          }else {
            $actual_masked -= $prevEnd - $start + 1;
          }
        }
      }
      $totalMasked += $actual_masked;
      $repeatsMasked++ if ( $actual_masked > 0 );

      my $origSeq = $maskSeqs{$seqID};
      substr($origSeq, $start-1, $len) = "N" x $len;
      $maskSeqs{$seqID} = $origSeq;
      
      if ( defined $INST ) {
        print $INST ""
            . $res->getScore . " "
            . $res->getQueryName() . " "
            . $res->getQueryStart() . " "
            . $res->getQueryEnd() . " "
      }

      $prevStart = $start;
      $prevEnd = $end if ( $end > $prevEnd);
      $prevID = $seqID;
    } # for
  } # if results > 0 
  undef $resultCollection;
  close $INST if ( defined $INST );

  # Write in input order so that the file is the same from run to run.
  open OUT, ">$fastaFile.masked";
  foreach my $seqID ( @maskSeqIDs ) {
      my $seq  = $maskSeqs{$seqID};
      print OUT ">$seqID\n";
      $seq =~ s/(.{50})/$1\n/g;
      print OUT "$seq\n";
  }
  close OUT;

  if ( $repeatsMasked == 0 ) {
    unlink( "$fastaFile.masked" );
    unlink( "$parameters{'instSeqFile'}" )
        if ( defined $parameters{'instSeqFile'}
             && -z $parameters{'instSeqFile'} );
  }
  unlink( grep { -e } $searchEngine->getSubjectArtifacts( $subject ) )
      unless ( $DEBUG );

  #print "    * Masked $repeatsMasked repeats totaling $totalMasked bp(s).\n";

  return ($repeatsMasked, $totalMasked);

}

#
# A helper function to open an input file, identify it as either a
#   "linup" =  Linup *.ali file
#   "msa-fasta" =  MSA file in FASTA format
#   "stockholm" =  Stockholm file
#   "malign" =  A *.malign file
#   "crossmatch" =  A crossmatch-like *.out file of one sequence vs many
#
# And return a fully populated MultAln object.
#
sub openAsMultAln{
  my $inputFile = shift;
  
  open IN, "<$inputFile" or die "openAsMultAln(): Could not open $inputFile for reading!\n";
  my $maxLines = 10000;
  my $fileType = "Unknown";
  my $foundFastaHdr = 0;
  my $foundConsensusHdr = 0;
  while ( <IN> )
  {
    next if (    /^\s*$/ 
              || /^(\W+).*Score:/ );
    last if ( $maxLines-- < 0 );
    if ( /^#\s+STOCKHOLM/ )
    {
      $fileType = "stockholm";
      last;
    }
    if ( /^\s*\d+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+\S+\s+\d+\s+\d+\s+\(\d+\)/ ||
         /Score:\s+\d+\s+Residues:/ )
    {
      # The second case is a check for a crossmatch file that has lots of binary
      # log entries preceeding the first alignment.
      $fileType = "crossmatch";
      last;
    }
    if ( /^consensus\s+\d+\s+\S+\s+\d+\s*$/ ) {
      $foundConsensusHdr = 1;
    }
    if ( $foundConsensusHdr && /^ref:/ ) {
      $fileType = "linup";
      last;
    }
    $foundFastaHdr = 1 if ( /^>\S+.*/ );
    if ( $foundFastaHdr && /^\s*[ACGTUMRWSYKVHDBNacgtumrwsykvhdbn\-\.]+\s*$/ )
    {
      $fileType = "msa-fasta";
      last;
    }
    if ( /^\s+\'alignCol\'\s+=>/ )
    {
      $fileType = "malign";
      last;
    }
      
  }
  close IN;
  
  if ( $fileType eq "Unknown" )
  {
    die "openAsMultAln(): Could not determine filetype for $inputFile.  Verify that\n"
        . "the file is either a cross_match, stockholm or an msa-fasta file.\n";
  }
  
  my $mAlign;
  my $seedAlign;
  if ( $fileType eq "crossmatch" )
  {
    my $resultCollection =
        CrossmatchSearchEngine::parseOutput( searchOutput => $inputFile );
  
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
      die "openAsMultAln(): Strange...this input file $inputFile appears not to be a multiple alignment!"
          if ( $staticQuery == 0 && $staticSubject == 0 );
      $queryID = $qID;
      $subjID  = $sID;
    }
    die "openAsMultAln(): Could not determine reference sequence in input file $inputFile.  This doesn't look like\n"
        . "a multiple alignment to one reference sequence!\n"
        if ( $staticQuery && $staticSubject );
  
    my $refInput = MultAln::Subject;
    $refInput = MultAln::Query if ( $staticQuery );
  
    $mAlign = MultAln->new(
                          referenceSeq              => "",
                          searchCollection          => $resultCollection,
                          searchCollectionReference => $refInput
                           );
  } elsif ( $fileType eq "stockholm" )
  {
    open my $IN, "<$inputFile" or die "openAsMultAln(): Could not open $inputFile for reading";
    $seedAlign = SeedAlignment->new();
    $seedAlign->read_stockholm( $IN );
    close $IN;
    $mAlign = MultAln->new( seedAlignment => $seedAlign );
  }elsif ( $fileType eq "msa-fasta" )
  { 
    my @seqs;
    my $seq;
    my $id;
    open my $IN, "<$inputFile" or die "openAsMultAln(): Could not open $inputFile for reading";
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
          warn "openAsMultAln(): WARN File contains a duplicate identifier \"$tmpID\".  A suffix of \"_$ver\"\n" .
               "                      will be appended to this occurence.\n";
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
    $mAlign = MultAln->new( sequences => \@seqs );
  }elsif ( $fileType eq "malign" ){
    $mAlign = MultAln->new();
    $mAlign = $mAlign->serializeIN( $inputFile );
  }elsif ( $fileType eq "linup" ) {
    open my $IN, "<$inputFile" or die "openAsMultAln(): Could not open $inputFile for reading";
    # Linup format
    my %seqHash = ();
    my $blockLen = 0;
    my $prefixWhitespace = 0;
    my $suffixWhitespace = 0;
    my $alignCols = 0;
    my $blockNumber = 0;
    my $newBlockHash = {};
    my $prevBlockHash;
    my %prevBlock = ();
    my %newBlock = ();
    my $nextLineIdx = 0;
    my $legacyFormat = 0;
    my $refSeq = "";
    my %uniqInStanza = ();
    while (<$IN>) {
      if ( /^consensus\s+\d+(\s+)(\S+)(\s+)/ ){
        $blockNumber++;
        $prefixWhitespace = length($1);
        $blockLen = length($2);
        $suffixWhitespace = length($3);
        $alignCols += $blockLen;
        $prevBlockHash = $newBlockHash;
        $newBlockHash = {};
        %uniqInStanza = ();
        next;
      }
      if ( /^ref:\S*\s+\d+\s{$prefixWhitespace}(.*)\s{$suffixWhitespace}\d+/ ) {
        # NOTE: This is working around a problem in which inexplicably some linup files
        # look like this:
        #
        #consensus                                           TACTGCTGACACGGAGAAAGTT 
        #ref:Tigger1#DNA/TcMar-Tigger                        TATTGCTGATATGGAGAAAGTT
        #fullTreeAnc111refChr24825_301065_302000_230_744_R   TACTGCTGACACGGAGAAAGTT
        #
        #consensus                                           NNNNNNNNNNNNNNNNNNNNNN 
        #ref:Tigger1#DNA/TcMar-Tigger                                    TAAGGAAAGA
        #fullTreeAnc111refChr21219_142236_144433_755_1394                TAAGGAAAGA
        #
        #consensus                                           TCCGTAACATAAAAGTGCANGG 
        #ref:Tigger1#DNA/TcMar-Tigger                        TCCATAACATAAAAGTGCAAGG
        #fullTreeAnc111refChr21219_142236_144433_755_1394    TCCGTAACATAAAAGTGCATGG
        #
        $refSeq .= $1;
        next;
      }
      if ( /^(\S+)\s+(\d+)\s{$prefixWhitespace}(.*)\s{$suffixWhitespace}(\d+)\s*$/ ) 
      {
        $legacyFormat = 1;
        # Legacy format
        #   In some cases a sequence may appear more than once in a MSA. For instance
        #   it could be due to a large deletion that caused an alignment break causing
        #   the same sequence (different regions ) to appear as independent fragments
        #   in the MSA.  In the traditional linup format this was difficult to parse
        #   as the identifier was not unique.  
        #
        #   Example 1 -- There are many legacy examples that have inconsistent
        #                indices.  In this example the alignment is in the forward
        #                direction yet the last stanza contains descending indices.
        #                Also of note...in the legacy output a stanza containing
        #                only gaps has a 1bp range as seen below.  The convention
        #                is that the start position refers the next nucleotide
        #                in the sequence ( e.g. start=1 means that the next
        #                nucleotide in the sequence will be at position 1). But
        #                what if no nucleotides exist in the stanza, only gap
        #                characters? I have chosen to have the start position point
        #                to the last nucleotide seen as in:
        #                "583 -- 583".  
        #              
        #   seq1    569 C---A--G--T--C--AC----C----AAA---T----A--C----A-    583
        #   seq1    584 --                                                  583
        #
        my $id = $1;
        my $start = $2;
        my $end = $4;
        my $orient = "+";
        $orient = "-" if ( $start > $end );
        my $seq = $3;
        my $lineID = "";
        # Attempt to join sequences between stanza's using ID + expected next start
        # position.
        if ( exists $prevBlockHash->{$id} ) {
          foreach my $rec ( @{$prevBlockHash->{$id}} ){
            # There are bugs in legacy Linup files
            # where the orientation appears to change if
            # the alignment ends on a block with only
            # "-" characters.
            #if ( $start == $rec->[0] && $orient eq $rec->[2] )
            if ( $start == $rec->[0] ){ 
              $lineID = $rec->[1];
            }
            if ( $start == ($rec->[0] - 1) ){
              die "RepeatUtil::openAsMultAln(): Legacy format with possible " . 
                  "inconsistent stanza joining [ id=$id, start=$start ].  Cannot parse.\n";
            }
          }
        }
        if ( exists $uniqInStanza{$id."_".$start} ) {
          die "RepeatUtil::openAsMultAln(): Legacy format with possible " . 
              "inconsistent stanza joining [ id=$id, start=$start ].  Cannot parse.\n";
        }
        $uniqInStanza{$id."_".$start}++;
        if ( ! exists $newBlockHash->{$id} ) {
          $newBlockHash->{$id} = [];
        }
        my $lid = $nextLineIdx;
        $lid = $lineID if ( $lineID ne "" );
        if ( $orient eq "+" ) {
          push @{$newBlockHash->{$id}}, [$end+1,$lid, $orient];
        }else {
          push @{$newBlockHash->{$id}}, [$end-1,$lid, $orient];
        }
 
        if ( $lineID eq "" ) {
          $lineID = $nextLineIdx;
          # Starting new sequence
          $seqHash{$lineID} = {};
          $seqHash{$lineID}->{'id'} = $id;
          $seqHash{$lineID}->{'start'} = $start;
          $seqHash{$lineID}->{'end'} = $end;
          $seqHash{$lineID}->{'orient'} = $orient;
          my $stanzaPadding = " "x($alignCols-$blockLen);
          $seqHash{$lineID}->{'seq'} = $stanzaPadding . $seq;
          $seqHash{$lineID}->{'lastblock'} = $blockNumber;
          $nextLineIdx++;
        }else {
          # Add to previous sequence
          $seqHash{$lineID}->{'end'} = $end; 
          if ( $end != $start ) {
            $seqHash{$lineID}->{'orient'} = $orient; 
          }
          $seqHash{$lineID}->{'seq'} .= $seq;
        }
      }

      # New LINUP format
      if ( /^(\S+)\s+(\d+)\s{$prefixWhitespace}(.*)\s{$suffixWhitespace}(\d+)\s+\[(\d+)\]$/ ) 
      {
        my $id = $1;
        my $start = $2;
        my $seq = $3;
        my $end = $4;
        my $orient = "+";
        $orient = "-" if ( $start > $end );
        my $lineID = $5;
        if ( exists $seqHash{$lineID} ) {
          # Append
          $seqHash{$lineID}->{'end'} = $end; 
          $seqHash{$lineID}->{'orient'} = $orient if ( $end != $start);
          $seqHash{$lineID}->{'seq'} .= $seq;
        }else {
          # New
          $seqHash{$lineID} = {};
          $seqHash{$lineID}->{'id'} = $id;
          $seqHash{$lineID}->{'start'} = $start;
          $seqHash{$lineID}->{'end'} = $end;
          $seqHash{$lineID}->{'orient'} = $orient; 
          my $stanzaPadding = " "x($alignCols-$blockLen);
          $seqHash{$lineID}->{'seq'} = $stanzaPadding . $seq;
          $seqHash{$lineID}->{'lastblock'} = $blockNumber;
        }
      }
    }

    if ( $legacyFormat ) {
      warn "WARNING: File $inputFile uses a legacy Linup format and may not import in a consistent fashion.\n";
    }
    my @seqs;
    foreach my $lineID ( sort {$a <=> $b} keys  %seqHash ) {
      my $id = $seqHash{$lineID}->{'id'};
      my $seq = $seqHash{$lineID}->{'seq'};
      my $start = $seqHash{$lineID}->{'start'};
      my $end = $seqHash{$lineID}->{'end'};
      my $orient = $seqHash{$lineID}->{'orient'};
      my $len = $end - $start + 1;
      # We now honor edge gaps if present
      #$seq =~ s/[-\s]/\./g;
      $seq =~ s/\-/\./g;
      push @seqs, [ $id, $seq, $start, $end ];
    }
    # Honoring edge gaps is for Linup format only ( for now and not optional )
    $mAlign = MultAln->new( reference => $refSeq, sequences => \@seqs, keepEdgeGaps => 1 );
  }else {
    die "openAsMultAln(): Support for $fileType is not complete yet ";
  }
  return($mAlign, $fileType);
}


1;





1;
