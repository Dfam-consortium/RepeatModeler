#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) TwoBitDB.pm
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      A read-only SeqDBI implementation over a UCSC 2bit file.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2026 Developed by
#* Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************

=head1 NAME

TwoBitDB - Read-only access to a UCSC 2bit sequence file

=head1 SYNOPSIS

use TwoBitDB;

  my $db = TwoBitDB->new( fileName     => "genome.2bit",
                          ucscToolsDir => "/usr/local/ucscTools" );

  my @ids     = $db->getIDs();
  my $len     = $db->getSeqLength( "gi|1" );
  my $nonAmb  = $db->getSubtLength( "gi|1" );
  my $seq     = $db->getSequence( "gi|1" );
  my $sub     = $db->getSubstr( "gi|1", 100, 50 );

  my $fh = $db->openRangeStream( [ { seqID => "gi|1", start => 1, end => 40000 },
                                   { seqID => "gi|2", start => 5, end => 900 } ] );
  while ( <$fh> ) { ... }
  close $fh;

=head1 DESCRIPTION

RepeatModeler keeps the genome it samples from in a UCSC 2bit file.  This
class reads such a file through the UCSC command line tools: sequence names
and lengths from twoBitInfo, ambiguous-base counts from twoBitInfo's N-block
output, and sequence from twoBitToFa.  It implements the read side of the
SeqDBI interface so that a caller written against FastaDB can take one of
these instead.

Bases come back upper case.  A 2bit file records lower case as a mask and
nothing else, and RepeatModeler treats case as meaningless in its samples,
so this class drops the mask rather than carry it into files that other
programs would then have to interpret.

The class writes every identifier a caller passes to a list file and hands
that to twoBitToFa through -seqList.  The -seq option cannot take a name
containing "/", and a RepeatMasker library is full of those.

Coordinates on the public methods follow SeqDBI: getSubstr() takes a
zero-based offset.  openRangeStream() takes one-based fully closed
coordinates because that is what RepeatModeler's sample source blocks carry.

=head1 SEE ALSO

SeqDBI, FastaDB

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

package TwoBitDB;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use File::Temp qw( tempfile );
use FindBin;
use lib $FindBin::RealBin;
use RepModelConfig;
use SeqDBI;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SeqDBI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "TwoBitDB";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##

=head2 new()

  Use: my $db = TwoBitDB->new( fileName     => "genome.2bit",
                               ucscToolsDir => "/usr/local/ucscTools",
                               [tempDir     => "/some/dir"] );

  ucscToolsDir defaults to the configured UCSCTOOLS_DIR.  tempDir is where
  the short-lived list files handed to twoBitToFa are written, and
  defaults to the system temporary directory.

=cut

##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  my $this = {};
  bless $this, $class;

  croak $CLASS . "::new(): Missing fileName parameter!\n"
      if ( !defined $nameValuePairs{'fileName'} );
  my $file = $nameValuePairs{'fileName'};
  croak $CLASS . "::new(): $file does not exist or is empty!\n"
      if ( !-s $file );
  croak $CLASS . "::new(): $file is not a 2bit file!\n"
      if ( !$class->isTwoBitFile( $file ) );

  my $toolsDir = $nameValuePairs{'ucscToolsDir'};
  $toolsDir = $RepModelConfig::configuration->{'UCSCTOOLS_DIR'}->{'value'}
      if ( !defined $toolsDir );
  foreach my $tool ( "twoBitInfo", "twoBitToFa" ) {
    croak $CLASS
        . "::new(): Cannot find $tool in $toolsDir.  Check the "
        . "UCSCTOOLS_DIR setting.\n"
        if ( !-x "$toolsDir/$tool" );
  }

  $this->{'fileName'}   = $file;
  $this->{'twoBitInfo'} = "$toolsDir/twoBitInfo";
  $this->{'twoBitToFa'} = "$toolsDir/twoBitToFa";
  $this->{'tempDir'}    = $nameValuePairs{'tempDir'} // File::Spec->tmpdir();

  $this->_loadIndex();

  return $this;
}

##-------------------------------------------------------------------------##

=head2 isTwoBitFile()

  Use: my $bool = TwoBitDB->isTwoBitFile( $path );

  True if $path starts with the 2bit signature, in either byte order.

=cut

##-------------------------------------------------------------------------##
sub isTwoBitFile {
  my $class = shift;
  my $path  = shift;

  return 0 if ( !defined $path || !-f $path );
  open my $fh, "<", $path or return 0;
  binmode $fh;
  my $magic;
  my $got = read( $fh, $magic, 4 );
  close $fh;
  return 0 if ( !defined $got || $got < 4 );

  return ( $magic eq pack( "N", 0x1A412743 )
           || $magic eq pack( "V", 0x1A412743 ) );
}

##-------------------------------------------------------------------------##
## Public Methods
##-------------------------------------------------------------------------##

sub getFileName {
  my $this = shift;

  return $this->{'fileName'};
}

sub getSeqCount {
  my $this = shift;

  return scalar( @{ $this->{'ids'} } );
}

sub getIDs {
  my $this = shift;

  return @{ $this->{'ids'} };
}

sub exists {
  my $this  = shift;
  my $seqID = shift;

  return ( exists $this->{'lengths'}->{$seqID} );
}

sub getSeqLength {
  my $this  = shift;
  my $seqID = shift;

  croak $CLASS . "::getSeqLength(): $seqID is not in " . $this->{'fileName'}
      if ( !exists $this->{'lengths'}->{$seqID} );

  return $this->{'lengths'}->{$seqID};
}

sub getSeqLengths {
  my $this = shift;

  return map { $this->{'lengths'}->{$_} } @{ $this->{'ids'} };
}

##-------------------------------------------------------------------------##

=head2 getXNLength()

  Use: my $count = $db->getXNLength( [$seqID] );

  The number of N bases in $seqID, or in the whole file when $seqID is
  omitted.  A 2bit file stores every ambiguity code as N, so this is also
  the count of every non-ACGT base the source held.

=cut

##-------------------------------------------------------------------------##
sub getXNLength {
  my $this  = shift;
  my $seqID = shift;

  $this->_loadNCounts();
  if ( defined $seqID ) {
    croak $CLASS . "::getXNLength(): $seqID is not in " . $this->{'fileName'}
        if ( !exists $this->{'lengths'}->{$seqID} );
    return $this->{'nCounts'}->{$seqID};
  }

  my $total = 0;
  $total += $this->{'nCounts'}->{$_} foreach ( @{ $this->{'ids'} } );
  return $total;
}

sub getXNLengths {
  my $this = shift;

  $this->_loadNCounts();
  return map { $this->{'nCounts'}->{$_} } @{ $this->{'ids'} };
}

sub getSubtLength {
  my $this  = shift;
  my $seqID = shift;

  if ( defined $seqID ) {
    return $this->getSeqLength( $seqID ) - $this->getXNLength( $seqID );
  }
  my $total = 0;
  $total += $this->{'lengths'}->{$_} foreach ( @{ $this->{'ids'} } );
  return $total - $this->getXNLength();
}

##-------------------------------------------------------------------------##

=head2 getNBlocks()

  Use: my @blocks = $db->getNBlocks( $seqID );

  The runs of N in $seqID as [ $start, $end ] pairs, zero-based and
  half-open, in ascending order.

=cut

##-------------------------------------------------------------------------##
sub getNBlocks {
  my $this  = shift;
  my $seqID = shift;

  $this->_loadNCounts();
  return @{ $this->{'nBlocks'}->{$seqID} // [] };
}

##-------------------------------------------------------------------------##

=head2 getSubtLengthInRange()

  Use: my $count = $db->getSubtLengthInRange( $seqID, $start, $end );

  Non-ambiguous bases in the one-based fully closed range $start..$end,
  computed from the N-block table without reading sequence.

=cut

##-------------------------------------------------------------------------##
sub getSubtLengthInRange {
  my $this  = shift;
  my $seqID = shift;
  my $start = shift;
  my $end   = shift;

  my $s0    = $start - 1;
  my $nBases = 0;
  foreach my $block ( $this->getNBlocks( $seqID ) ) {
    last if ( $block->[ 0 ] >= $end );
    next if ( $block->[ 1 ] <= $s0 );
    my $a = $block->[ 0 ] > $s0 ? $block->[ 0 ] : $s0;
    my $b = $block->[ 1 ] < $end ? $block->[ 1 ] : $end;
    $nBases += $b - $a;
  }

  return ( $end - $start + 1 ) - $nBases;
}

sub getDescription {
  my $this = shift;

  # A 2bit file holds names only.
  return "";
}

sub getDescriptors {
  my $this = shift;

  return map { "" } @{ $this->{'ids'} };
}

sub getSequence {
  my $this  = shift;
  my $seqID = shift;

  croak $CLASS . "::getSequence(): $seqID is not in " . $this->{'fileName'}
      if ( !exists $this->{'lengths'}->{$seqID} );

  return $this->_extractOne( $seqID );
}

##-------------------------------------------------------------------------##

=head2 getSubstr()

  Use: my $sequence = $db->getSubstr( $seqID, $offset, [$length] );

  $offset is zero-based.  Without $length everything from $offset to the
  end of the sequence is returned.

=cut

##-------------------------------------------------------------------------##
sub getSubstr {
  my $this   = shift;
  my $seqID  = shift;
  my $offset = shift;
  my $length = shift;

  my $seqLen = $this->getSeqLength( $seqID );
  $length = $seqLen - $offset if ( !defined $length );
  return "" if ( $length <= 0 || $offset >= $seqLen );
  my $end = $offset + $length;
  $end = $seqLen if ( $end > $seqLen );

  return $this->_extractOne( "$seqID:$offset-$end" );
}

##-------------------------------------------------------------------------##

=head2 openRangeStream()

  Use: my $fh = $db->openRangeStream( \@ranges );

  Extract many ranges in one twoBitToFa run and return a filehandle
  reading FASTA records in the order requested.  Each range is a hash
  with seqID, start and end, where start and end are one-based and fully
  closed.  A range without start and end selects the whole sequence.

  twoBitToFa names each record "seqID:start-end" in its own zero-based
  half-open convention.  Callers that care about the names should
  re-header by position, which is what sampleFromDB() does.

=cut

##-------------------------------------------------------------------------##
sub openRangeStream {
  my $this   = shift;
  my $ranges = shift;

  my ( $lfh, $listFile ) =
      tempfile( "twobitdb-XXXXXX", DIR => $this->{'tempDir'}, SUFFIX => ".lst" );
  foreach my $range ( @{$ranges} ) {
    croak $CLASS
        . "::openRangeStream(): "
        . $range->{'seqID'}
        . " is not in "
        . $this->{'fileName'} . "\n"
        if ( !exists $this->{'lengths'}->{ $range->{'seqID'} } );
    if ( defined $range->{'start'} && defined $range->{'end'} ) {
      print $lfh $range->{'seqID'} . ":"
          . ( $range->{'start'} - 1 ) . "-"
          . $range->{'end'} . "\n";
    }
    else {
      print $lfh $range->{'seqID'} . "\n";
    }
  }
  close $lfh;

  return $this->_openExtraction( "-seqList=$listFile", $listFile );
}

##-------------------------------------------------------------------------##

=head2 openFastaStream()

  Use: my $fh = $db->openFastaStream();

  A filehandle reading the whole file as FASTA, in file order.

=cut

##-------------------------------------------------------------------------##
sub openFastaStream {
  my $this = shift;

  return $this->_openExtraction( "" );
}

##-------------------------------------------------------------------------##
## Unsupported SeqDBI write methods
##-------------------------------------------------------------------------##

sub addSequence {
  croak $CLASS . " is read-only.\n";
}

sub removeSequence {
  croak $CLASS . " is read-only.\n";
}

sub setSubstr {
  croak $CLASS . " is read-only.\n";
}

sub compact {
  return;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

sub _loadIndex {
  my $this = shift;

  my $cmd = $this->{'twoBitInfo'} . " " . $this->{'fileName'} . " stdout";
  open my $fh, "$cmd |"
      or croak $CLASS . "::_loadIndex(): Could not run $cmd: $!\n";
  my @ids;
  my %lengths;
  while ( <$fh> ) {
    next unless ( /^(\S+)\t(\d+)/ );
    push @ids, $1;
    $lengths{$1} = $2;
  }
  close $fh;
  croak $CLASS
      . "::_loadIndex(): twoBitInfo reported no sequences in "
      . $this->{'fileName'}
      . ".  The file may be damaged.\n"
      if ( $? != 0 || !@ids );

  $this->{'ids'}     = \@ids;
  $this->{'lengths'} = \%lengths;
}

sub _loadNCounts {
  my $this = shift;

  return if ( exists $this->{'nCounts'} );

  my %counts = map { $_ => 0 } @{ $this->{'ids'} };
  my %blocks;
  my $cmd = $this->{'twoBitInfo'} . " -nBed " . $this->{'fileName'} . " stdout";
  open my $fh, "$cmd |"
      or croak $CLASS . "::_loadNCounts(): Could not run $cmd: $!\n";
  while ( <$fh> ) {
    next unless ( /^(\S+)\t(\d+)\t(\d+)/ );
    $counts{$1} += $3 - $2;
    push @{ $blocks{$1} }, [ $2, $3 ];
  }
  close $fh;
  croak $CLASS . "::_loadNCounts(): $cmd failed.\n" if ( $? != 0 );

  $this->{'nCounts'} = \%counts;
  $this->{'nBlocks'} = \%blocks;
}

# Run twoBitToFa with the given selection option and return a filehandle
# on its output.  $cleanupFile, if given, is removed when the handle is
# closed.
sub _openExtraction {
  my $this        = shift;
  my $selection   = shift;
  my $cleanupFile = shift;

  my $cmd = $this->{'twoBitToFa'} . " -noMask $selection "
      . $this->{'fileName'} . " stdout";
  my $pid = open( my $fh, "-|", $cmd );
  croak $CLASS . "::_openExtraction(): Could not run $cmd: $!\n"
      if ( !$pid );

  # Tie the temporary file's lifetime to the handle's.
  $this->{'pending'}->{ fileno( $fh ) } = $cleanupFile
      if ( defined $cleanupFile );

  return $fh;
}

# Read one record from twoBitToFa and return its bases.
sub _extractOne {
  my $this = shift;
  my $spec = shift;

  my ( $lfh, $listFile ) =
      tempfile( "twobitdb-XXXXXX", DIR => $this->{'tempDir'}, SUFFIX => ".lst" );
  print $lfh "$spec\n";
  close $lfh;

  my $fh  = $this->_openExtraction( "-seqList=$listFile", $listFile );
  my $seq = "";
  while ( <$fh> ) {
    next if ( /^>/ );
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  $this->_closeExtraction( $fh );

  return $seq;
}

sub _closeExtraction {
  my $this = shift;
  my $fh   = shift;

  my $file = delete $this->{'pending'}->{ fileno( $fh ) };
  close $fh;
  unlink( $file ) if ( defined $file && -e $file );
}

##-------------------------------------------------------------------------##

=head2 closeStream()

  Use: $db->closeStream( $fh );

  Close a handle from openRangeStream() or openFastaStream() and remove
  the list file behind it.  A plain close() works too; the list file is
  then removed when the object is destroyed.

=cut

##-------------------------------------------------------------------------##
sub closeStream {
  my $this = shift;
  my $fh   = shift;

  $this->_closeExtraction( $fh );
}

sub DESTROY {
  my $this = shift;

  foreach my $file ( values %{ $this->{'pending'} // {} } ) {
    unlink( $file ) if ( defined $file && -e $file );
  }
}

1;
