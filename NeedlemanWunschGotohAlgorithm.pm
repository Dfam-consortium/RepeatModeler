#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) NeedlemanWunschGotohAlgorithm.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A perl implementation of the Needleman/Wunsch/Gotoh
##      global alignment algorithm.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless( {}
#     , 'NeedlemanWunschGotohAlgorithm' );
#
###############################################################################
# ChangeLog
#
#     $Log: NeedlemanWunschGotohAlgorithm.pm,v $
#     Revision 1.22  2017/04/05 00:03:31  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#
#

=head1 NAME

NeedlemanWunschGotohAlgorithm

=head1 SYNOPSIS

use NeedlemanWunschGotohAlgorithm

Usage: 


    my $querySeq = "AACACTCCATGGATGGGGGGTTAGGATGCGGCATTAT";
    my $subjSeq = "AACACCCCATGGATGGTTAGGATGCGGCATTAT";
    my $searchResult = NeedlemanWunschGotohAlgorithm::search(
                                  querySeq => $querySeq, 
                                  subjectSeq => $subjSeq,
                                  matrixFile => "/Matrices/linupmatrix", 
                                  insOpenPenalty => -12, 
                                  insExtPenalty => -2,
                                  delOpenPenalty => -12,  
                                  delExtPenalty => -2 );

    print "Result:\n" . 
          $searchResult->toStringFormatted( SearchResult::AlignWithQuerySeq ) . 
          "\n";                                                                                


=head1 DESCRIPTION

=head1 INSTANCE METHODS

=cut 

package NeedlemanWunschGotohAlgorithm;
use strict;
use Carp;
use SearchResult;
use Data::Dumper;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "NeedlemanWunschGotohAlgorithm";

my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 search()

  Use: my $searchResult = search( 
                        subjectSeq => "AAAAAAA",
                        querySeq => "GAAAGT",
                        matrix => SequencSimilarityMatrix::new,
                        insOpenPenalty => -25,
                        insExtPenalty => -5,
                        delOpenPenalty => -25,
                        delExtPenalty => -5 
                                );

  A simple implementation of the Needleman-Wunsch-Gotoh Algorithm
  as described in Biological Sequence Analysis ( R. Durbin, S. Eddy,
  A. Krogh, and G. Mitchison ). This is a global sequence alignment
  algorithm with support for affine gap penalties, and alignment
  output.

=cut

##-------------------------------------------------------------------------##
sub search {
  my %parameters = @_;

  croak $CLASS
      . "::needlemanWunschGotohSearch(): Missing subjectSeq "
      . "parameter!\n"
      unless ( defined $parameters{'subjectSeq'} );
  my $subjSeq = $parameters{'subjectSeq'};

  croak $CLASS
      . "::needlemanWunschGotohSearch(): Missing querySeq "
      . "parameter!\n"
      unless ( defined $parameters{'querySeq'} );
  my $querySeq = $parameters{'querySeq'};

  croak $CLASS
      . "::needlemanWunschGotohSearch(): Missing matrix "
      . "parameter!\n"
      unless ( defined $parameters{'matrix'} );
  my $matrixRef = $parameters{'matrix'};

  croak $CLASS
      . "::needlemanWunschGotohSearch(): Missing insOpenPenalty "
      . "parameter!\n"
      unless ( defined $parameters{'insOpenPenalty'} );
  my $insOpenPenalty = $parameters{'insOpenPenalty'};

  croak $CLASS
      . "::needlemanWunschGotohSearch(): Missing insExtPenalty "
      . "parameter!\n"
      unless ( defined $parameters{'insExtPenalty'} );
  my $insExtPenalty = $parameters{'insExtPenalty'};

  croak $CLASS
      . "::needlemanWunschGotohSearch(): Missing delOpenPenalty "
      . "parameter!\n"
      unless ( defined $parameters{'delOpenPenalty'} );
  my $delOpenPenalty = $parameters{'delOpenPenalty'};

  croak $CLASS
      . "::needlemanWunschGotohSearch(): Missing delExtPenalty "
      . "parameter!\n"
      unless ( defined $parameters{'delExtPenalty'} );
  my $delExtPenalty = $parameters{'delExtPenalty'};

  #print "Aligning $querySeq to $subjSeq\n";

  my $score          = 0;
  my $queryAlignment = "";
  my $subjAlignment  = "";

  my @dynMatrix = ();
  $dynMatrix[ 0 ][ 0 ] = {
                           'ins' => 0,
                           'del' => 0,
                           'sub' => 0
  };
  my $subFromIns;
  my $subFromDel;
  my $subFromSub;
  my $insFromIns;
  my $insFromSub;
  my $delFromDel;
  my $delFromSub;
  my $matrixValue;
  my $penalty;

  my $i;
  my $j;

  my @queryArray = split( //, uc( $querySeq ) );
  my @subjArray  = split( //, uc( $subjSeq ) );

  for ( $i = 0 ; $i <= length( $subjSeq ) ; $i++ ) {
    for ( $j = 0 ; $j <= length( $querySeq ) ; $j++ ) {
      next if ( $i == 0 && $j == 0 );
      if ( $i == 0 ) {
        $penalty = $delOpenPenalty + ( ( $j - 1 ) * $delExtPenalty );
        $dynMatrix[ 0 ][ $j ] = {
                                  'sub' => -10000,
                                  'ins' => -10000,
                                  'del' => $penalty
        };
        next;
      }
      if ( $j == 0 ) {
        $penalty = $insOpenPenalty + ( ( $i - 1 ) * $insExtPenalty );
        $dynMatrix[ $i ][ 0 ] = {
                                  'sub' => -10000,
                                  'ins' => $penalty,
                                  'del' => -10000
        };
        next;
      }

# TODO Check with Arian on this
#$matrixValue =
#    $matrixRef->{'matrixHash'}->{ $queryArray[ $j - 1 ], $subjArray[ $i - 1 ] };
      $matrixValue =
          $matrixRef->{'matrixHash'}
          ->{ $subjArray[ $i - 1 ], $queryArray[ $j - 1 ] };

      # Calculate sub value
      $subFromSub = $dynMatrix[ $i - 1 ][ $j - 1 ]->{'sub'};
      $subFromIns = $dynMatrix[ $i - 1 ][ $j - 1 ]->{'ins'};
      $subFromDel = $dynMatrix[ $i - 1 ][ $j - 1 ]->{'del'};

      if ( $subFromIns >= $subFromSub && $subFromIns >= $subFromDel ) {

        # Insertion
        $dynMatrix[ $i ][ $j ]->{'sub'}          = $subFromIns + $matrixValue;
        $dynMatrix[ $i ][ $j ]->{'tracebackSub'} = "ins";
      }
      elsif ( $subFromSub > $subFromDel && $subFromSub > $subFromIns ) {

        # Substition
        $dynMatrix[ $i ][ $j ]->{'sub'}          = $subFromSub + $matrixValue;
        $dynMatrix[ $i ][ $j ]->{'tracebackSub'} = "sub";
      }
      else {

        # Deletion
        $dynMatrix[ $i ][ $j ]->{'sub'}          = $subFromDel + $matrixValue;
        $dynMatrix[ $i ][ $j ]->{'tracebackSub'} = "del";
      }

      # Calculate ins value
      $insFromIns = $dynMatrix[ $i - 1 ][ $j ]->{'ins'} + $insExtPenalty;
      $insFromSub = $dynMatrix[ $i - 1 ][ $j ]->{'sub'} + $insOpenPenalty;
      if ( $insFromSub >= $insFromIns ) {

        # Substition
        $dynMatrix[ $i ][ $j ]->{'ins'}          = $insFromSub;
        $dynMatrix[ $i ][ $j ]->{'tracebackIns'} = "sub";
      }
      else {

        # Insertion
        $dynMatrix[ $i ][ $j ]->{'ins'}          = $insFromIns;
        $dynMatrix[ $i ][ $j ]->{'tracebackIns'} = "ins";
      }

      # Calculate del value
      $delFromDel = $dynMatrix[ $i ][ $j - 1 ]->{'del'} + $delExtPenalty;
      $delFromSub = $dynMatrix[ $i ][ $j - 1 ]->{'sub'} + $delOpenPenalty;
      if ( $delFromSub >= $delFromDel ) {

        # Substition
        $dynMatrix[ $i ][ $j ]->{'del'}          = $delFromSub;
        $dynMatrix[ $i ][ $j ]->{'tracebackDel'} = "sub";
      }
      else {

        # Deletion
        $dynMatrix[ $i ][ $j ]->{'del'}          = $delFromDel;
        $dynMatrix[ $i ][ $j ]->{'tracebackDel'} = "del";
      }
    }
  }

  if ( $DEBUG ) {
    print "" . ( $i * $j ) . " cells calculated\n";
    print " i = $i j = $j  first = " . $#dynMatrix . "\n";
    print "dynMatrix = \n$subjSeq\n";
    for ( my $l = 0 ; $l <= length( $subjSeq ) ; $l++ ) {

      for ( my $k = 0 ; $k <= length( $querySeq ) ; $k++ ) {
        print "\t\t" . $dynMatrix[ $l ][ $k ]->{'sub'};
      }
      print "\n";

      for ( my $k = 0 ; $k <= length( $querySeq ) ; $k++ ) {
        print "\t\t" . $dynMatrix[ $l ][ $k ]->{'ins'};
      }
      print "\n";

      for ( my $k = 0 ; $k <= length( $querySeq ) ; $k++ ) {
        print "\t\t" . $dynMatrix[ $l ][ $k ]->{'del'};
      }
      print "\n";

      print "\n\n";
    }
  }

  # Set pointers to the last cell of the matrix
  $i--;
  $j--;

  my $lastEmission = "";
  if (    $dynMatrix[ $i ][ $j ]->{'sub'} >= $dynMatrix[ $i ][ $j ]->{'del'}
       && $dynMatrix[ $i ][ $j ]->{'sub'} >= $dynMatrix[ $i ][ $j ]->{'ins'} )
  {
    $score        = $dynMatrix[ $i ][ $j ]->{'sub'};
    $lastEmission = "sub";
  }
  elsif (    $dynMatrix[ $i ][ $j ]->{'ins'} > $dynMatrix[ $i ][ $j ]->{'sub'}
          && $dynMatrix[ $i ][ $j ]->{'ins'} > $dynMatrix[ $i ][ $j ]->{'del'} )
  {
    $score        = $dynMatrix[ $i ][ $j ]->{'ins'};
    $lastEmission = "ins";
  }
  else {
    $score        = $dynMatrix[ $i ][ $j ]->{'del'};
    $lastEmission = "del";
  }

  while ( !( $i == 0 || $j == 0 ) ) {

    if ( $lastEmission eq "del" ) {
      $lastEmission = $dynMatrix[ $i ][ $j ]->{'tracebackDel'};

      # Deletion
      $queryAlignment .= substr( $querySeq, $j - 1, 1 );
      $subjAlignment .= "-";
      $j--;
    }
    elsif ( $lastEmission eq "sub" ) {
      $lastEmission = $dynMatrix[ $i ][ $j ]->{'tracebackSub'};
      $queryAlignment .= substr( $querySeq, $j - 1, 1 );
      $subjAlignment  .= substr( $subjSeq,  $i - 1, 1 );
      $i--;
      $j--;

    }
    elsif ( $lastEmission eq "ins" ) {
      $lastEmission = $dynMatrix[ $i ][ $j ]->{'tracebackIns'};

      # Insertion
      $queryAlignment .= "-";
      $subjAlignment .= substr( $subjSeq, $i - 1, 1 );
      $i--;

    }
  }

  if ( 0 ) {
    while ( !( $i == 0 || $j == 0 ) ) {

      if (  $dynMatrix[ $i ][ $j ]->{'del'} >= $dynMatrix[ $i ][ $j ]->{'sub'}
         && $dynMatrix[ $i ][ $j ]->{'del'} >= $dynMatrix[ $i ][ $j ]->{'ins'} )
      {

        # Deletion
        $queryAlignment .= substr( $querySeq, $j - 1, 1 );
        $subjAlignment .= "-";
        $j--;
      }
      elsif ( $dynMatrix[ $i ][ $j ]->{'sub'} >= $dynMatrix[ $i ][ $j ]->{'del'}
         && $dynMatrix[ $i ][ $j ]->{'sub'} >= $dynMatrix[ $i ][ $j ]->{'ins'} )
      {
        $queryAlignment .= substr( $querySeq, $j - 1, 1 );
        $subjAlignment  .= substr( $subjSeq,  $i - 1, 1 );
        $i--;
        $j--;
      }
      else {

        # Insertion
        $queryAlignment .= "-";
        $subjAlignment .= substr( $subjSeq, $i - 1, 1 );
        $i--;
      }

    }

  }
  if ( $i != 0 ) {
    $subjAlignment .= reverse( substr( $subjSeq, 0, $i ) );
    $queryAlignment .= "-" x ( $i );
  }
  if ( $j != 0 ) {
    $queryAlignment .= reverse( substr( $querySeq, 0, $j ) );
    $subjAlignment .= "-" x ( $j );
  }

  $queryAlignment = reverse( $queryAlignment );
  $subjAlignment  = reverse( $subjAlignment );

  my $result = SearchResult->new(
                                  queryName      => "query",
                                  subjName       => "subject",
                                  pctInsert      => 0.0,
                                  pctDelete      => 0.0,
                                  subjStart      => 1,
                                  subjEnd        => length( $subjSeq ),
                                  queryStart     => 1,
                                  queryEnd       => length( $querySeq ),
                                  score          => $score,
                                  pctDiverge     => 0.0,
                                  subjRemaining  => 0,
                                  queryRemaining => 0,
                                  orientation    => "",
                                  queryString    => $queryAlignment,
                                  subjString     => $subjAlignment
  );

  return $result;

}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

## TODO: Merge this and the other many matrix objects into one!
sub readMatrix {
  my ( $matrixFileName ) = @_;

  open MATRIX, "<$matrixFileName" || die "Can't open $matrixFileName!\n";

  my $dataRow  = 0;
  my %matrix   = ();
  my @alphabet = ();
  while ( <MATRIX> ) {
    $_ = uc( $_ );
    chomp;
    next if ( /FREQS\s+(.*)/ );
    next if ( /^#/ );

    if ( /^\s*[A-Z]\s+[A-Z]\s+[A-Z]\s+[A-Z]\s+/ ) {
      @alphabet = split;
    }
    elsif ( @alphabet
            && /^\s*(\W\s+)*[\d-]+\s+[\d-]+\s+[\d-]+\s+[\d-]+\s+/ )
    {
      my @rowValues = split;
      if ( $1 ne "" ) {
        shift @rowValues;
      }
      for ( my $i = 0 ; $i < $#rowValues ; $i++ ) {
        $matrix{ $alphabet[ $dataRow ] }->{ $alphabet[ $i ] } =
            $rowValues[ $i ];
      }
      $dataRow++;
    }
  }
  close MATRIX;

  return \%matrix;
}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = toString([$this]);
##
##      $this         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
sub toString {
  my $this = shift;
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##	  $filename	: A filename to be created
##
##  Returns
##
##	Uses the Data::Dumper module to save out the data
##	structure as a text file.  This text file can be
##	read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $this     = shift;
  my $fileName = shift;

  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

##-------------------------------------------------------------------------##
## Use: my serializeIN( $filename );
##
##	$filename	: A filename containing a serialized object
##
##  Returns
##
##	Uses the Data::Dumper module to read in data
##	from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $this         = shift;
  my $fileName     = shift;
  my $fileContents = "";
  my $oldSep       = $/;
  undef $/;
  my $in;
  open $in, "$fileName";
  $fileContents = <$in>;
  $/            = $oldSep;
  close $in;
  return eval( $fileContents );
}

1;

