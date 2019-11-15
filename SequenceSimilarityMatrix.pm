#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) SequenceSimilarityMatrix.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A generic similarity matrix reader
##
#******************************************************************************
#*  Copyright (C) 2001 by Robert Hubley, Institute for Systems Biology.       *
#*  All rights reserved.                                                      *
#*                                                                            *
#*  This software is part of a set of software utilities for processing       *
#*  of genetic sequence data.  It should not be redistributed or used for     *
#*  any commercial purpose, including commercially funded sequencing,         *
#*  without written permission from the author and the Institute for          *
#*  Systems Biology.                                                          *
#*                                                                            *
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology be     *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#*                                                                            *
#******************************************************************************
#bless( {
#                 'monoAlphabetHash' => {
#                                         'N' => 1,
#                                         .
#                                         'G' => 1,
#                                       },
#                 'monoSubTotals' => {
#                                      'G' => '106347',
#                                      .
#                                      'T' => '183590'
#                                    },
#                 'monoSubCounts' => {
#                                      'GN' => '1',
#                                      .
#                                      'CT' => '26405'
#                                    },
#                 'monoSubProb' => {
#                                    'TT' => '0.843362928264067',
#                                    .
#                                    'TC' => '0.074328667138733'
#                                  },
#                 'monoSubHash' => {
#                                    'TT' => '0.91730283955087',
#                                    .
#                                    'CT' => '0.139273651004344'
#                                  },
#                 'monoSubDivergence' => '0.139085568540752',
#                 'monoSubScale' => '0.480595582861076',
#                 'monoRateHash' => {
#                                     'TT' => '-0.18872880068887',
#                                     .
#                                     'CT' => '0.336182624927569'
#                                   },
#                 'monoEquilHash' => {
#                                      'TT' => '0.373615931039398',
#                                      .
#                                      'CT' => '0.373615931036906'
#                                    },
#                 'diSubCounts' => {
#                                    'GCGG' => '459',
#                                    .
#                                    'TTCC' => '331'
#                                  },
#                 'diEquilHash' => {
#                                    'CCGT' => '0.0474291179601828',
#                                    .
#                                    'GCGC' => '0.0189083125022889'
#                                  },
#                 'diSubProb' => {
#                                  'GGCA' => '0.0153687079945481',
#                                  .
#                                  'CTCG' => '0.00460998830800067'
#                                },
#                 'triSubProb' => {
#                                   'TAGCGC' => '0.00046029919447641',
#                                   .
#                                   'AATTCT' => '0.000962721292185912'
#                                 },
#                 'trippleSubTotals' => {
#                                         'AAA' => '30242',
#                                         .
#                                         'TTT' => '30242'
#                                       },
#                 'triSubScale' => '0.359724715125594',
#                 'triEquilHash' => {
#                                     'ACGGGA' => '0.00587170136071818',
#                                     .
#                                     'AATTCT' => '0.0189427945712668'
#                                   },
#                 'triSubDivergence' => '0.138205411835254',
#                 'triSubHash' => {
#                                   'ACGGGA' => '-0.00124954148882776',
#                                   .
#                                   'AATTCT' => '0.000185848349802029'
#                                 },
#                 'trippleAlphabetHash' => {
#                                            'AAA' => 1,
#                                            .
#                                            'TTT' => 1
#                                          },
#                 'triSubPosCounts' => {
#                                        'AAA' => {
#                                                   '7491' => {
#                                                               'G' => '3',
#                                                               'A' => '204',
#                                                               'C' => '2',
#                                                               'T' => '6'
#                                                             },
#                                                   '1176' => {
#                                                               'G' => '10',
#                                                               'A' => '237',
#                                                               'C' => '12',
#                                                               'T' => '13'
#                                                             }
#                                        'NAA' => {
#                                                   '3097' => {
#                                                               'G' => '19',
#                                                               'A' => '130',
#                                                               'C' => '4',
#                                                               'T' => '18'
#                                                             }
#                                                 }
#                                      },
#                 'triRateHash' => {
#                                    'ACGGGA' => '-0.0123740187206802',
#                                    .
#                                    'AATTCT' => '0.000230782000486845'
#                                  },
#                 'trippleSubCounts' => {
#                                         'TAGCGC' => '1',
#                                         .
#                                         'AATTCT' => '17'
#                                       }
#               }, 'SequenceSimilarityMatrix' )
#
###############################################################################
#
# ChangeLog
#
#   $Log: SequenceSimilarityMatrix.pm,v $
#   Revision 1.22  2017/04/05 00:03:31  rhubley
#   Cleanup before a distribution
#
#
###############################################################################
# To Do:
#
#

=head1 NAME

SequenceSimilarityMatrix

=head1 SYNOPSIS

use SequenceSimilarityMatrix;

Usage:

    my $matrix = SequenceSimilarityMatrix->new();

=head1 DESCRIPTION

A class for storing a generic sequence similarity
matrix.

=head1 SEE ALSO

=over 4

  RepeatModeler

=back

=head1 COPYRIGHT

Copyright 2004 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut

package SequenceSimilarityMatrix;
use strict;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "SequenceSimilarityMatrix";

##-------------------------------------------------------------------------##
##
##
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;
  my $this  = {};
  bless $this, $class;
  return $this;
}

##---------------------------------------------------------------------##

=head2 parseFromFile()

  Use: $obj->parseFromFile( "/usr/lib/matrices/blosom64" );

  Read a matrix file into the object.

=cut

##---------------------------------------------------------------------##
sub parseFromFile {
  my $object   = shift;
  my $fileName = shift
      || croak $CLASS. "::parseFromFile() needs a filename!";

  open MATRIX, $fileName
      or croak "Unable to open matrix file: $fileName : $!";

  my %freqHash  = ();
  my @rowValues = ();
  my $row       = 0;
  my $i         = 0;
  while ( <MATRIX> ) {
    $_ = uc( $_ );
    chomp;
    if ( /FREQS\s+(.*)/ ) {
      %freqHash = split " ", $1;
    }
    if ( /^\s*[A-Z]\s+[A-Z]\s+[A-Z]\s+[A-Z]\s+/ ) {
      chomp;
      s/\s//g;
      $object->{alphabet} = $_;
      $object->{alphabetArray} = [ split( '', $_ ) ];
    }
    if ( /^\s*([a-zA-Z]\s+)?[\d-]+\s+[\d-]+\s+[\d-]+\s+[\d-]+\s+/ ) {
      @rowValues = split;
      shift @rowValues if ( $1 ne "" );
      my $col = 0;
      foreach $a ( @{ $object->{alphabetArray} } ) {
        $object->{matrixHash}{ $object->{alphabetArray}[ $row ], $a } =
            $rowValues[ $col++ ];
      }
      $object->{matrix}[ $row++ ] = [ @rowValues ];
    }
  }
  close MATRIX;

  if ( scalar( keys %freqHash ) == 0 ) {
    $freqHash{"A"} = 0.25;
    $freqHash{"C"} = 0.25;
    $freqHash{"G"} = 0.25;
    $freqHash{"T"} = 0.25;
  }

  for ( $i = 0 ; $i < length( $object->{alphabet} ) ; $i++ ) {
    $object->{freqArray}[ $i ] =
        $freqHash{ substr( $object->{alphabet}, $i, 1 ) } || 0;
  }

}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: toString([$object]);
##
##      $object         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
## TODO: This could be moved into a generic object class.
##       or placed in perlhelpers module.
sub toString {
  my $object = shift;
  my $data_dumper = new Data::Dumper( [ $object ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##        $filename     : A filename to be created
##
##  Returns
##
##      Uses the Data::Dumper module to save out the data
##      structure as a text file.  This text file can be
##      read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $object   = shift;
  my $fileName = shift;

  my $data_dumper = new Data::Dumper( [ $object ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

##-------------------------------------------------------------------------##
## Use: my serializeIN( $filename );
##
##      $filename       : A filename containing a serialized object
##
##  Returns
##
##      Uses the Data::Dumper module to read in data
##      from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $object       = shift;
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

##-------------------------------------------------------------------------##

=head2 addMatrices()

  Use:   my  addMatrices( $matrixObject,
                          $addCounts,
                          $addRates,
                          $addSubs );


     $matrixObject  : Add the counts (and alphabet) from
                      $matrixObject to this one.
     $addCounts     : Boolean.  Add mono/tri counts.
     $addRates      : Boolean.  Add mono/tri rates.
     $addSubs       : Boolean.  Add mono/tri substitution values.


  Modifies the object by adding the data from another
  matrix object to this one.  NOTE: It also adds to the
  alphabet if necessary.

=cut

##-------------------------------------------------------------------------##
sub addMatrices {
  my $object    = shift;
  my $matrix    = shift;
  my $addCounts = shift;
  my $addRates  = shift;
  my $addSubs   = shift;

  # Mono nucleotide counts
  foreach my $cBase ( keys( %{ $matrix->{monoAlphabetHash} } ) ) {

    # Add to alphabets if necessary
    $object->{monoAlphabetHash}{$cBase} = 1;

    # Add matrix counts to object counts
    foreach my $sBase ( keys( %{ $matrix->{monoAlphabetHash} } ) ) {
      $object->{monoSubCounts}{ $cBase . $sBase } +=
          $matrix->{monoSubCounts}{ $cBase . $sBase }
          if ( $addCounts && $matrix->{monoSubCounts}{ $cBase . $sBase } > 0 );
      $object->{monoRateHash}{ $cBase . $sBase } +=
          $matrix->{monoRateHash}{ $cBase . $sBase }
          if ( $addRates && $matrix->{monoRateHash}{ $cBase . $sBase } > 0 );
      $object->{monoSubHash}{ $cBase . $sBase } +=
          $matrix->{monoSubHash}{ $cBase . $sBase }
          if ( $addSubs && $matrix->{monoSubHash}{ $cBase . $sBase } > 0 );
    }
  }

  # Tri nucleotide counts
  foreach my $cBase ( keys( %{ $matrix->{trippleAlphabetHash} } ) ) {

    # Add to alphabets if necessary
    $object->{trippleAlphabetHash}{$cBase} = 1;

    # Add matrix counts to object counts
    foreach my $sBase ( keys( %{ $matrix->{trippleAlphabetHash} } ) ) {
      $object->{trippleSubCounts}{ $cBase . $sBase } +=
          $matrix->{trippleSubCounts}{ $cBase . $sBase }
          if (    $addCounts
               && $matrix->{trippleSubCounts}{ $cBase . $sBase } > 0 );
      $object->{triRateHash}{ $cBase . $sBase } +=
          $matrix->{triRateHash}{ $cBase . $sBase }
          if ( $addRates && $matrix->{triRateHash}{ $cBase . $sBase } > 0 );
      $object->{triSubHash}{ $cBase . $sBase } +=
          $matrix->{triSubHash}{ $cBase . $sBase }
          if ( $addSubs && $matrix->{triSubHash}{ $cBase . $sBase } > 0 );
    }
  }
}

##-------------------------------------------------------------------------##

=head2 div()

  Use:   my  div($value,
                 $counts,
                 $rates,
                 $subs );

       $value      : The value to divide the matrix by
       $counts     : Boolean.  divide mono/tri counts by value.
       $rates      : Boolean.  divide mono/tri rates by value.
       $subs       : Boolean.  divide mono/tri substitution by value.

  Modifies the object by dividing the matrix by a constant
  value.

=cut

##-------------------------------------------------------------------------##
sub divMatrices {
  my $object = shift;
  my $value  = shift;
  my $counts = shift;
  my $rates  = shift;
  my $subs   = shift;

  # Mono nucleotide matrices
  foreach my $cBase ( keys( %{ $object->{monoAlphabetHash} } ) ) {
    foreach my $sBase ( keys( %{ $object->{monoAlphabetHash} } ) ) {
      $object->{monoSubCounts}{ $cBase . $sBase } /= $value
          if ( $counts && $object->{monoSubCounts}{ $cBase . $sBase } > 0 );
      $object->{monoRateHash}{ $cBase . $sBase } /= $value
          if ( $rates && $object->{monoRateHash}{ $cBase . $sBase } > 0 );
      $object->{monoSubHash}{ $cBase . $sBase } /= $value
          if ( $subs && $object->{monoSubHash}{ $cBase . $sBase } > 0 );
    }
  }

  # Tri nucleotide matrices
  foreach my $cBase ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
    foreach my $sBase ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
      $object->{trippleSubCounts}{ $cBase . $sBase } /= $value
          if ( $counts && $object->{trippleSubCounts}{ $cBase . $sBase } > 0 );
      $object->{triRateHash}{ $cBase . $sBase } /= $value
          if ( $rates && $object->{triRateHash}{ $cBase . $sBase } > 0 );
      $object->{triSubHash}{ $cBase . $sBase } /= $value
          if ( $subs && $object->{triSubHash}{ $cBase . $sBase } > 0 );
    }
  }
}

##-------------------------------------------------------------------------##
## Use:   my  $divergence = _matrixDivergence( \@matrixArrayRef,
##                                             $alphabetRef,
##                                             [$monoFreqHashRef] );
##
##          \@matrixArrayRef      : A square probability matrix
##          $AlphabetArrayRef     : Array header
##          $monoFreqHashRef      : A hash of mono nucleotides and
##                                  their frequencies within the
##                                  background for which this matrix
##                                  is tuned.
##
##      Returns
##
##          Performs one of two types of divergence calculations
##          based on the presence of mono-nucleotide background
##          frequencies.
##
##          If frequencies are not present it calculates the
##          divergence as the inverse of the average diagonal
##          probability.
##
##          If the frequences are present it calculates the
##          divergence as the inverse of the sum of weighted
##          probabilities.  ie.
##
##               1 - ( p(A)*p(A->A) + p(C) * p(C->C)...)
##
##-------------------------------------------------------------------------##
sub _matrixDivergence {
  my $matrixArrayRef   = shift;
  my $alphabetArrayRef = shift;
  my $monoFreqHashRef  = shift;

  my $diagTotal = 0;
  for ( my $i = 0 ; $i <= $#$matrixArrayRef ; $i++ ) {
    if ( defined $monoFreqHashRef ) {
      $diagTotal +=
          $matrixArrayRef->[ $i ][ $i ] *
          $monoFreqHashRef->{ $alphabetArrayRef->[ $i ] };
    }
    else {
      $diagTotal += $matrixArrayRef->[ $i ][ $i ];
    }
  }
  if ( defined $monoFreqHashRef ) {
    return ( 1 - ( $diagTotal ) );
  }
  else {
    return ( 1 - ( $diagTotal / ( $#$matrixArrayRef + 1 ) ) );
  }
}

##-------------------------------------------------------------------------##
## Use:   my  $divergence = _triMatrixDivergence( $matrixArrayRef
##                                                $alphabetArrayRef,
##                                                [$noCGSites],
##                                                [$trippletFreqHashRef]);
##
##        $matrixArrayRef        : A square probability matrix
##        $alphabetArrayRef      : The row/column identifiers
##        $noCGSites             : Optional parameter to ignore
##                                 sites containing CGs.
##        $trippletFreqHashRef   : Optional parameter to use
##                                 tripplet backgound frequencies to
##                                 weight the calculation.
##
## Returns
##
##          Performs one of several types of divergence calculations
##          based on the presence of tri-nucleotide background
##          frequencies (trippletFreqHashRef), and or the flag
##          $noCGSites.
##
##          The basic divergence calculation (when neither $noCGSites
##          or $trippletFreqHashRef are set/defined) is defined as
##          the inverse average of all loose-perfect probabilities.
##          Ie.
##            1 - Avg(p(AAA->NAN) + p(AAT->NAN) + p(AAC->NAN) +....)
##
##          The $noCGSites option simply removed all p(NCG->NCN) and
##          p(CGN->NGN) terms from the previous calculation.
##
##          If $trippletFreqHashRef is defined a very different
##          formula is used. NOTE: That the $noCGSites option
##          can also be used with this formula but may not
##          make much sense:
##
##            1 - (p(AAA) * p(AAA->NAN) + p(AAT) * p(AAT->NAN) + ....)
##
##-------------------------------------------------------------------------##
sub _triMatrixDivergence {
  my $matrixArrayRef      = shift;
  my $alphabetArrayRef    = shift;
  my $noCGSites           = shift;
  my $trippletFreqHashRef = shift;

  my $diagTotal    = 0;
  my $diagElements = 0;
  for ( my $i = 0 ; $i <= $#$matrixArrayRef ; $i++ ) {
    next if ( $noCGSites && index( $alphabetArrayRef->[ $i ], "CG" ) >= 0 );
    $diagElements++;
    for ( my $j = 0 ; $j <= $#$matrixArrayRef ; $j++ ) {
      if (
           substr( $alphabetArrayRef->[ $i ], 1, 1 ) eq
           substr( $alphabetArrayRef->[ $j ], 1, 1 ) )
      {
        if ( defined $trippletFreqHashRef ) {
          $diagTotal +=
              ( $matrixArrayRef->[ $i ][ $j ] *
                $trippletFreqHashRef->{ $alphabetArrayRef->[ $i ] } );
        }
        else {
          $diagTotal += $matrixArrayRef->[ $i ][ $j ];
        }
      }
    }
  }
  if ( defined $trippletFreqHashRef ) {
    return ( 1 - ( $diagTotal ) );
  }
  else {
    return ( 1 - ( $diagTotal / $diagElements ) );
  }
}

sub hashToArray {
  my $object               = shift;
  my $matrixHashRef        = shift;
  my $alphabetArrayRefOrig = shift;
  my ( @matrixArray, $alphabetArrayRef ) =
      _hashToArray( $matrixHashRef, $alphabetArrayRefOrig );
  return ( \@matrixArray, $alphabetArrayRef );
}

##-------------------------------------------------------------------------##
## Use:   my  \@matrixArray = _hashToArray( \%matrixHashRef,
##					     [$alphabetArrayRef] );
##
##              \%matrixHashRef  : A hash which contains array information
##				   stored as two concatenated equal length
##				   keys.  ie.
##
##					A B       %matrixHash{AA} = 0
##				      A 0 1       %matrixHash{AB} = 1
##				      B 1 0       %matrixHash{BA} = 1
##                                                %matrixHash{BB} = 0
##
##		$alphabetArrayRef : Optional.  Instead of using the
##				    the keys from the hash to determine
##				    the matrix dimensions...use this
##				    array as key pairs in the hash.
##
##      Returns
##
##		A PERL array form of the hash matrix data.
##
##-------------------------------------------------------------------------##
sub _hashToArray {
  my $matrixHashRef    = shift;
  my $alphabetArrayRef = shift;

  # Calculate the unique set of key pairs for the
  # hash if we are not given this information already
  unless ( defined $alphabetArrayRef ) {
    my %seen = ();
    $alphabetArrayRef = [
      grep    { !$seen{$_}++ }
          map {
        (
          substr( $_, 0,                length( $_ ) / 2 ),
          substr( $_, length( $_ ) / 2, length( $_ ) / 2 )
            )
          }
          keys( %{$matrixHashRef} )
    ];
  }

  # Build the array.
  my $i           = 0;
  my $j           = 0;
  my @matrixArray = ();
  foreach my $conSeq ( @{$alphabetArrayRef} ) {
    foreach my $subSeq ( @{$alphabetArrayRef} ) {
      $matrixArray[ $i ][ $j++ ] = $matrixHashRef->{ $conSeq . $subSeq };
    }
    $j = 0;
    $i++;
  }

  return ( \@matrixArray, $alphabetArrayRef );
}

##-------------------------------------------------------------------------##
## Use:   my  \%matrixHash = _arrayToHash( \@matrixArrayRef,
##					   \@alphabetArrayRef );
##
##              \@matrixArrayef    : An square array
##		\@alphabetArrayRef : The col/row identifiers to be
##				     used as keys in the hash.
##
##				     ie.
##
##					A B       %matrixHash{AA} = 0
##				      A 0 1       %matrixHash{AB} = 1
##				      B 1 0       %matrixHash{BA} = 1
##                                                %matrixHash{BB} = 0
##
##      Returns
##
##		A PERL hash form of the array matrix data.
##
##-------------------------------------------------------------------------##
sub _arrayToHash {
  my $matrixArrayRef   = shift;
  my $alphabetArrayRef = shift;

  # Build the hash.
  my $i          = 0;
  my $j          = 0;
  my %matrixHash = ();
  foreach my $conSeq ( @{$alphabetArrayRef} ) {
    foreach my $subSeq ( @{$alphabetArrayRef} ) {
      $matrixHash{ $conSeq . $subSeq } = $matrixArrayRef->[ $i ][ $j++ ];
    }
    $j = 0;
    $i++;
  }
  return ( \%matrixHash );
}

##-------------------------------------------------------------------------##

=head2 matrixToString()

  TODO: DOCUMENT!

=cut

##-------------------------------------------------------------------------##
sub matrixToString {
  my $object     = shift;
  my $name       = shift;
  my $compressed = shift;
  my $excludeNs  = shift;
  return _matrixHashToString( $object->{$name}, $compressed, $excludeNs );
}

##-------------------------------------------------------------------------##

=head2 triStandardDeviation()

  TODO: DOCUMENT!

=cut

##-------------------------------------------------------------------------##
sub triStandardDeviation {
  my $object = shift;
  my $outStr = "";

  my %seen = ();
  my @alphabet = grep { !$seen{$_}++; }
      map { ( substr( $_, 0, 3 ), substr( $_, 3 ) ); }
      keys( %{ $object->{trippleSubCounts} } );
  %seen = ();
  my @sortedMidAlphabet =
      ( sort grep { !$seen{$_}++; }
        ( map { ( substr( $_, 1, 1 ) ); } @alphabet ) );
  my @sortedAlphabet = (
    sort {
             substr( $a, 1, 1 ) cmp substr( $b, 1, 1 )
          || substr( $a, 0, 1 ) cmp substr( $b, 0, 1 )
          || substr( $a, 2, 1 ) cmp substr( $b, 2, 1 )
        } @alphabet
  );
  $outStr .= "Tripplet Position Standard Deviation\n";
  $outStr .= "====================================\n";
  $outStr .= "\t";
  foreach my $tripplet ( ( "A", "C", "G", "T" ) ) {
    $outStr .= "$tripplet\t";
  }
  $outStr .= "\n";
  foreach my $tripplet ( @sortedAlphabet ) {
    next if ( ( $tripplet =~ /N/ ) );
    my %triBaseTotals        = ();
    my $numOfPositionsInSeq  = 0;
    my %meanBaseCount        = ();
    my %sumOfSquaredDistance = ();
    foreach my $pos ( keys( %{ $object->{triSubPosCounts}{$tripplet} } ) ) {
      $numOfPositionsInSeq++;
      foreach
          my $sub ( keys( %{ $object->{triSubPosCounts}{$tripplet}{$pos} } ) )
      {
        next if ( $sub eq "N" );
        $triBaseTotals{$sub} +=
            $object->{triSubPosCounts}{$tripplet}{$pos}{$sub};
        $meanBaseCount{$sub} = $triBaseTotals{$sub} / $numOfPositionsInSeq++;
      }
    }
    foreach my $pos ( keys( %{ $object->{triSubPosCounts}{$tripplet} } ) ) {
      foreach
          my $sub ( keys( %{ $object->{triSubPosCounts}{$tripplet}{$pos} } ) )
      {
        next if ( $sub eq "N" );
        $sumOfSquaredDistance{$sub} +=
            ( $object->{triSubPosCounts}{$tripplet}{$pos}{$sub} -
              $meanBaseCount{$sub} )**2;
      }
    }
    $outStr .= "$tripplet\t";
    foreach my $sub ( sort keys( %sumOfSquaredDistance ) ) {
      $object->{sdTriBaseCounts}{$tripplet}{$sub} +=
          sqrt(
           ( 1 / ( $numOfPositionsInSeq - 1 ) ) * $sumOfSquaredDistance{$sub} );
      $outStr .=
          sprintf( "%0.6g", $object->{sdTriBaseCounts}{$tripplet}{$sub} )
          . "\t";
    }
    $outStr .= "\n";
  }

  return ( $outStr );

}

##-------------------------------------------------------------------------##
##
##
##-------------------------------------------------------------------------##
sub _matrixHashToString {
  my $matrixHashRef = shift;
  my $compressed    = shift;
  my $excludeNs     = shift;
  my $outStr        = "";

  #
  #  Sort the matrix col/row headings (alphabet)
  #
  my @sortedAlphabet    = ();
  my @sortedMidAlphabet = ();
  my @alphabet          = ();
  my %seen              = ();
  my $keyLen            = length( ( keys( %{$matrixHashRef} ) )[ 0 ] );
  if ( $excludeNs ) {
    @alphabet = grep { !$seen{$_}++ && !/N/; }
        map { ( substr( $_, 0, $keyLen / 2 ), substr( $_, $keyLen / 2 ) ); }
        keys( %{$matrixHashRef} );
  }
  else {
    @alphabet = grep { !$seen{$_}++; }
        map { ( substr( $_, 0, $keyLen / 2 ), substr( $_, $keyLen / 2 ) ); }
        keys( %{$matrixHashRef} );
  }

  if ( $keyLen == 6 ) {
    my %seen = ();
    @sortedMidAlphabet =
        ( sort grep { !$seen{$_}++; }
          ( map { ( substr( $_, 1, 1 ) ); } @alphabet ) );
    @sortedAlphabet = (
      sort {
               substr( $a, 1, 1 ) cmp substr( $b, 1, 1 )
            || substr( $a, 0, 1 ) cmp substr( $b, 0, 1 )
            || substr( $a, 2, 1 ) cmp substr( $b, 2, 1 )
          } @alphabet
    );
  }
  else {
    @sortedAlphabet = ( sort( @alphabet ) );
  }

  #
  # Form the report string
  #
  $outStr .= "\t";
  my $colIndex = "";
  my %colTot   = ();

  # Generate the column headings
  if ( $compressed ) {
    foreach $colIndex ( @sortedMidAlphabet ) {
      $outStr .= "$colIndex\t";
    }
  }
  else {
    foreach $colIndex ( @sortedAlphabet ) {
      $outStr .= "$colIndex\t";
    }
  }
  $outStr .= "\n";

  # Generate the rows
  my $rowIndex = "";
  foreach $rowIndex ( @sortedAlphabet ) {
    $outStr .= "$rowIndex\t";
    foreach $colIndex ( @sortedAlphabet ) {
      if ( $compressed ) {
        $colTot{ substr( $colIndex, 1, 1 ) } +=
            $matrixHashRef->{ $rowIndex . $colIndex };
      }
      else {
        $outStr .= ""
            . sprintf( "%.6g", $matrixHashRef->{ $rowIndex . $colIndex } )
            . "\t";
      }
    }
    if ( $compressed ) {
      map { $outStr .= sprintf( "%.6g", $colTot{$_} ) . "\t"; }
          @sortedMidAlphabet;
      %colTot = undef;
    }
    $outStr .= "\n";
  }
  $outStr .= "\n\n";

  return ( $outStr );
}

##-------------------------------------------------------------------------##

=head2 calculateLLRMatrices()

  TODO: DOCUMENT!

=cut

##-------------------------------------------------------------------------##
sub calculateLLRMatrices {
  my $object       = shift;
  my $backgroundGC = shift;
  my $scale        = shift;

  my $outStr = "FREQS A "
      . ( ( 1 - $backgroundGC ) / 2 ) . " C "
      . ( $backgroundGC / 2 ) . " G "
      . ( $backgroundGC / 2 ) . " T "
      . ( ( 1 - $backgroundGC ) / 2 ) . "\n";
  $outStr .= "A\tR\tG\tC\tY\tT\tK\tM\tS\tW\tN\tX\n";

  #
  # Create a IUB Matrix and a background vector
  #
  my @cBaseArray = ();
  my @sBaseArray = ();
  foreach my $cBase ( "A", "R", "G", "C", "Y", "T", "K", "M", "S", "W" ) {
    if ( $cBase eq "A" || $cBase eq "T" || $cBase eq "G" || $cBase eq "C" ) {
      @cBaseArray = ( $cBase );
    }
    elsif ( $cBase eq "R" ) {

      # A or G (Purines)
      @cBaseArray = ( "A", "G" );
    }
    elsif ( $cBase eq "Y" ) {

      # C or T (Pyrimidines)
      @cBaseArray = ( "C", "T" );
    }
    elsif ( $cBase eq "K" ) {

      # G or T (Keto)
      @cBaseArray = ( "G", "T" );
    }
    elsif ( $cBase eq "M" ) {

      # A or C (Amino)
      @cBaseArray = ( "A", "C" );
    }
    elsif ( $cBase eq "S" ) {

      # G or C (Strong 3H bonds)
      @cBaseArray = ( "G", "C" );
    }
    elsif ( $cBase eq "W" ) {

      # A or T (Weak 2H bonds)
      @cBaseArray = ( "A", "T" );
    }
    elsif ( $cBase eq "N" ) {

      # Any Base
    }
    elsif ( $cBase eq "X" ) {

      # Mask Character
    }
    foreach my $sBase ( "A", "R", "G", "C", "Y", "T", "K", "M", "S", "W" ) {
      if ( $sBase eq "A" || $sBase eq "T" || $sBase eq "G" || $sBase eq "C" ) {
        @sBaseArray = ( $sBase );
      }
      elsif ( $sBase eq "R" ) {

        # A or G (Purines)
        @sBaseArray = ( "A", "G" );
      }
      elsif ( $sBase eq "Y" ) {

        # C or T (Pyrimidines)
        @sBaseArray = ( "C", "T" );
      }
      elsif ( $sBase eq "K" ) {

        # G or T (Keto)
        @sBaseArray = ( "G", "T" );
      }
      elsif ( $sBase eq "M" ) {

        # A or C (Amino)
        @sBaseArray = ( "A", "C" );
      }
      elsif ( $sBase eq "S" ) {

        # G or C (Strong 3H bonds)
        @sBaseArray = ( "G", "C" );
      }
      elsif ( $sBase eq "W" ) {

        # A or T (Weak 2H bonds)
        @sBaseArray = ( "A", "T" );
      }
      elsif ( $sBase eq "N" ) {

        # Any Base
      }
      elsif ( $sBase eq "X" ) {

        # Mask Character
      }
      my $foregroundProb = 0;
      my $backgroundProb = 0;
      foreach my $firstBase ( @cBaseArray ) {
        foreach my $secondBase ( @sBaseArray ) {
          $foregroundProb += $object->{monoSubHash}{ $firstBase . $secondBase };
          if ( $secondBase eq "A" || $secondBase eq "T" ) {
            $backgroundProb += ( 1 - $backgroundGC ) / 2;
          }
          else {
            $backgroundProb += $backgroundGC / 2;
          }
        }
      }
      $outStr .= sprintf( "%1.0f",
                          ( log( $foregroundProb ) - log( $backgroundProb ) ) *
                              $scale )
          . "\t";
    }
    $outStr .= "\n";
  }
  return ( $outStr );
}

##-------------------------------------------------------------------------##

=head2 calculateProbMatrix

  TODO: DOCUMENT

=cut

##-------------------------------------------------------------------------##
#
# TODO: Currently excludes only N's in the matrix.....hmmmm
# TODO: Also...do we really need to store the ..SubTotals in the
#       datastructure?
#
sub calculateProbMatrix {
  my $object = shift;

  #
  #
  #
  foreach my $cBase ( keys( %{ $object->{monoAlphabetHash} } ) ) {
    next if ( $cBase eq "N" );
    foreach my $sBase ( keys( %{ $object->{monoAlphabetHash} } ) ) {
      next if ( $sBase eq "N" );
      $object->{monoSubTotals}{$cBase} +=
          $object->{monoSubCounts}{ $cBase . $sBase };
    }
    foreach my $sBase ( keys( %{ $object->{monoAlphabetHash} } ) ) {
      next if ( $sBase eq "N" );
      if ( $object->{monoSubTotals}{$cBase} > 0 ) {
        $object->{monoSubProb}{ $cBase . $sBase } =
            $object->{monoSubCounts}{ $cBase . $sBase } /
            $object->{monoSubTotals}{$cBase};
      }
    }
  }

  #
  # Tripplet prob using pseudo counts (laplace's method)
  #
  # TODO: Remove trippleSubTotals from the data structure
  #       since we don't need it.
  my $totalTripplets = 0;
  foreach my $cTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
    next if ( index( $cTripplet, "N" ) >= 0 );
    foreach my $sTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
      next if ( index( $sTripplet, "N" ) >= 0 );
      $object->{trippleSubTotals}{$cTripplet} +=
          $object->{trippleSubCounts}{ $cTripplet . $sTripplet } + 1;
    }
    $totalTripplets += $object->{trippleSubTotals}{$cTripplet};
    foreach my $sTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
      next if ( index( $sTripplet, "N" ) >= 0 );
      if ( $object->{trippleSubTotals}{$cTripplet} > 0 ) {
        $object->{triSubProb}{ $cTripplet . $sTripplet } =
            ( $object->{trippleSubCounts}{ $cTripplet . $sTripplet } + 1 ) /
            $object->{trippleSubTotals}{$cTripplet};
      }
    }
  }

}

##-------------------------------------------------------------------------##

=head2 printMatrixReport()

  TODO: DOCUMENT!

=cut

##-------------------------------------------------------------------------##
sub printMatrixReport {
  my $object          = shift;
  my $compressTriData = shift;
  my $outStr          = "";

  $outStr .= "Mono-nucleotide count matrix:\n";
  $outStr .= "=============================\n";
  $outStr .= _matrixHashToString( $object->{monoSubCounts} );

  $outStr .= "Mono-nucleotide probability matrix:\n";
  $outStr .= "===================================\n";
  $outStr .= _matrixHashToString( $object->{monoSubProb} );

  $outStr .= "Mono-nucleotide rate matrix:\n";
  $outStr .= "============================\n";
  $outStr .= _matrixHashToString( $object->{monoRateHash} );

  $outStr .= "Mono-nucleotide substitution matrix:\n";
  $outStr .= "====================================\n";
  $outStr .= "\tscale = $object->{monoSubScale}\n";
  $outStr .= "\tdivergence = $object->{monoSubDivergence}\n\n";
  $outStr .= _matrixHashToString( $object->{monoSubHash} );

  $outStr .= "Tri-nucleotide count matrix:\n";
  $outStr .= "============================\n";
  $outStr .=
      _matrixHashToString( $object->{trippleSubCounts}, $compressTriData );

  $outStr .= "Tri-nucleotide probability matrix:\n";
  $outStr .= "==================================\n";
  $outStr .= _matrixHashToString( $object->{triSubProb}, $compressTriData );

  $outStr .= "Tri-nucleotide rate matrix:\n";
  $outStr .= "===========================\n";
  $outStr .= _matrixHashToString( $object->{triRateHash}, $compressTriData );

  $outStr .= "Tri-nucleotide substitution matrix:\n";
  $outStr .= "===================================\n";
  $outStr .= "\tscale = $object->{triSubScale}\n";
  $outStr .= "\tdivergence = $object->{triSubDivergence}\n\n";
  $outStr .= _matrixHashToString( $object->{triSubHash}, $compressTriData );

  return ( $outStr );

}

##-------------------------------------------------------------------------##

=head2 generateFlankingBaseGraphData()

 Use: my $outStr = generateFlankingBaseGraphData( $subBase )

       $subBase :    The base to consider flanking
                     effects on.

  A string with a tab delimited probability matrix.
  The matrix is of the form:

                      A   C   T
                 AGA  .01 .03 .04
                 AGT  .02 .01 .01
                 ...


=cut

##-------------------------------------------------------------------------##
sub generateFlankingBaseGraphData {
  my $object              = shift;
  my $subBase             = shift;    # Base to consider substitions on
  my $fixed               = shift;    # Fixed or loose flanking bases
  my $toString            = shift;    # Create a string or struct output
  my $singleStrand        = shift;    # Single stranded counts?
  my $optTrippletArrayRef = shift;    # Set order for data

  my @trippletAlphabet = keys( %{ $object->{trippleAlphabetHash} } );
  if ( defined $optTrippletArrayRef ) {
    @trippletAlphabet = @$optTrippletArrayRef;
  }

  #
  # combine data to create a *X*->Y mapping
  #   Where X is the base of interest (parameter $subBase)
  #   and Y is one of the three other bases it could be
  #   substituted for.
  #
  my %outArr = ();
  my $outStr = "\t";
  foreach my $sBase ( "A", "G", "C", "T" ) {
    next if ( $sBase eq $subBase );
    $outStr .= $sBase . "\t";
  }
  $outStr .= "\t\tTotal Datapoints\n";
  foreach my $cTripplet ( @trippletAlphabet ) {
    next
        if ( index( $cTripplet, "N" ) >= 0
             || substr( $cTripplet, 1, 1 ) ne $subBase );
    my %sTrippletCounts = ();
    my $sTrippletTotal  = 0;
    foreach my $sTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
      next if ( index( $sTripplet, "N" ) >= 0 );
      next
          if (
               $fixed
               && (   substr( $sTripplet, 0, 1 ) ne substr( $cTripplet, 0, 1 )
                   || substr( $sTripplet, 2, 1 ) ne substr( $cTripplet, 2, 1 ) )
          );
      my $sBase = substr( $sTripplet, 1, 1 );
      if ( $singleStrand ) {
        $sTrippletCounts{$sBase} +=
            $object->{trippleSubSingleCounts}{ $cTripplet . $sTripplet };
        $sTrippletTotal +=
            $object->{trippleSubSingleCounts}{ $cTripplet . $sTripplet };
      }
      else {
        $sTrippletCounts{$sBase} +=
            $object->{trippleSubCounts}{ $cTripplet . $sTripplet };
        $sTrippletTotal +=
            $object->{trippleSubCounts}{ $cTripplet . $sTripplet };
      }
    }
    $outStr .= "$cTripplet\t";
    foreach my $sBase ( "A", "G", "C", "T" ) {
      next if ( $sBase eq $subBase );
      $outStr .= ""
          . ( ( $sTrippletCounts{$sBase} + 1 ) / ( $sTrippletTotal + 1 ) )
          . "\t";
      push(
            @{ $outArr{$cTripplet} },
            ( ( $sTrippletCounts{$sBase} + 1 ) / ( $sTrippletTotal + 1 ) )
      );

    }
    $outStr .= "\t\t$sTrippletTotal";
    $outStr .= "\n";
  }

  # Return the results
  if ( $toString ) {
    return $outStr;
  }
  else {
    return \%outArr;
  }
}

##-------------------------------------------------------------------------##

=head2 generateFlankingBaseConsistencyGraphData()

 Use:   my $outStr = generateFlankingBaseConsistencyGraphData(
                                                 $subBase,
                                                 $toString,
                                                 $optTrippletArrayRef)

  NOTE: This is different from generateFlankingBaseGraphData
        in that it uses the 14% substitution matrices in order
        to compare more than one element against each other.

=cut

##-------------------------------------------------------------------------##
sub generateFlankingBaseConsistencyGraphData {
  my $object              = shift;
  my $subBase             = shift;    # Base to consider substitions on
  my $toString            = shift;    # Create a string or struct output
  my $optTrippletArrayRef = shift;    # Set order for data

  my @trippletAlphabet = keys( %{ $object->{trippleAlphabetHash} } );
  if ( defined $optTrippletArrayRef ) {
    @trippletAlphabet = @$optTrippletArrayRef;
  }

  #
  # combine data to create a *X*->Y mapping
  #   Where X is the base of interest (parameter $subBase)
  #   and Y is one of the three other bases it could be
  #   substituted for.
  #
  my %outArr = ();
  my $outStr = "\t";
  foreach my $sBase ( "A", "G", "C", "T" ) {
    next if ( $sBase eq $subBase );
    $outStr .= $sBase . "\t";
  }
  $outStr .= "\n";
  foreach my $cTripplet ( @trippletAlphabet ) {
    next
        if ( index( $cTripplet, "N" ) >= 0
             || substr( $cTripplet, 1, 1 ) ne $subBase );
    my %sTrippletCounts = ();
    foreach my $sTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
      next if ( index( $sTripplet, "N" ) >= 0 );
      my $sBase = substr( $sTripplet, 1, 1 );
      $sTrippletCounts{$sBase} +=
          $object->{triSubHash}{ $cTripplet . $sTripplet };
    }
    $outStr .= "$cTripplet\t";
    foreach my $sBase ( "A", "G", "C", "T" ) {
      next if ( $sBase eq $subBase );
      $outStr .= "" . $sTrippletCounts{$sBase} . "\t";
      push( @{ $outArr{$cTripplet} }, $sTrippletCounts{$sBase} );
    }
    $outStr .= "\n";
  }

  # Return the results
  if ( $toString ) {
    return $outStr;
  }
  else {
    return \%outArr;
  }
}

##-------------------------------------------------------------------------##

=head2 generateTriPosGraphData()

  Use:   my $outStr = generateTriPosGraphData( $tripplet )

=cut

##-------------------------------------------------------------------------##
sub generateTriPosGraphData {
  my $object   = shift;
  my $tripplet = shift;

  my $outStr = "ACGT\n";
  foreach my $triPos ( keys( %{ $object->{'triSubPosCounts'}{$tripplet} } ) ) {
    $outStr .= $triPos . "\t";
    $outStr .= $object->{'triSubPosCounts'}{$tripplet}{$triPos}{'A'} . "\t";
    $outStr .= $object->{'triSubPosCounts'}{$tripplet}{$triPos}{'C'} . "\t";
    $outStr .= $object->{'triSubPosCounts'}{$tripplet}{$triPos}{'G'} . "\t";
    $outStr .= $object->{'triSubPosCounts'}{$tripplet}{$triPos}{'T'} . "\t";
    $outStr .= "\n";
  }
  return $outStr;

}

##-------------------------------------------------------------------------##

=head2 generateMonoFromTri()

  Use:   my $outStr = generateMonoFromTri( $excludeCpG )

  Mono substitution frequencies with or without
  CpG consensus sites.  In order to use flanking
  consensus base information the Mono substitution
  frequencies need to be generated from the tri
  nucleotide data.

=cut

##-------------------------------------------------------------------------##
sub generateMonoFromTriDeleteMe {
  my $object     = shift;
  my $excludeCpG = shift;

  my %monoSubProb   = ();
  my %monoSubCounts = ();
  my %monoSubTotals = ();
  foreach my $cTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
    next if ( index( $cTripplet, "N" ) >= 0 );
    next if ( $excludeCpG == 1 && index( $cTripplet, "CG" ) >= 0 );
    foreach my $sTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
      next if ( index( $sTripplet, "N" ) >= 0 );
      $monoSubCounts{ substr( $cTripplet, 1, 1 )
            . substr( $sTripplet, 1, 1 ) } +=
          $object->{trippleSubCounts}{ $cTripplet . $sTripplet };
      $monoSubTotals{ substr( $cTripplet, 1, 1 ) } +=
          $object->{trippleSubCounts}{ $cTripplet . $sTripplet };
    }
  }
  foreach my $sBase ( "GA", "AG", "GT", "AC", "AT", "GC" ) {
    my $cplBase = $sBase;
    $cplBase =~ tr/ACGTRYKMSWBDHV/TGCAYRMKWSVHDB/;
    $monoSubProb{ $sBase . " "
          . $cplBase } = (
       (
         ( $monoSubCounts{$sBase} / $monoSubTotals{ substr( $sBase, 0, 1 ) } ) +
             (
           $monoSubCounts{$cplBase} / $monoSubTotals{ substr( $cplBase, 0, 1 ) }
             )
       ) / 2
          );
  }

  return ( \%monoSubProb );
}

##-------------------------------------------------------------------------##

=head2 generateMonoFromTri()

  Use:   my $outStr = generateMonoFromTri( $excludeCpG )

  Mono substitution frequencies with or without
  CpG consensus sites.  In order to use flanking
  consensus base information the Mono substitution
  frequencies need to be generated from the tri
  nucleotide data.

=cut

##-------------------------------------------------------------------------##
sub generateMonoFromTri {
  my $object           = shift;
  my $excludeCpG       = shift;
  my $thisTrippletOnly = shift;
  my $triFreqRef       = shift;

  my %monoSubProb   = ();
  my %monoSubCounts = ();
  my %monoSubTotals = ();
  foreach my $cTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
    next if ( index( $cTripplet, "N" ) >= 0 );
    next if ( $excludeCpG == 1 && index( $cTripplet, "CG" ) >= 0 );
    next if ( $thisTrippletOnly ne "" && $cTripplet ne $thisTrippletOnly );
    foreach my $sTripplet ( keys( %{ $object->{trippleAlphabetHash} } ) ) {
      next if ( index( $sTripplet, "N" ) >= 0 );
      if ( defined $triFreqRef ) {
        $monoSubCounts{ substr( $cTripplet, 1, 1 )
              . substr( $sTripplet, 1, 1 ) } +=
            ( $object->{trippleSubCounts}{ $cTripplet . $sTripplet } *
              $triFreqRef->{$cTripplet} );
        $monoSubTotals{ substr( $cTripplet, 1, 1 ) } +=
            ( $object->{trippleSubCounts}{ $cTripplet . $sTripplet } *
              $triFreqRef->{$cTripplet} );
      }
      else {
        $monoSubCounts{ substr( $cTripplet, 1, 1 )
              . substr( $sTripplet, 1, 1 ) } +=
            $object->{trippleSubCounts}{ $cTripplet . $sTripplet };
        $monoSubTotals{ substr( $cTripplet, 1, 1 ) } +=
            $object->{trippleSubCounts}{ $cTripplet . $sTripplet };
      }
    }
  }
  foreach my $sBase ( "GA", "AG", "GT", "AC", "AT", "GC" ) {
    my $cplBase = $sBase;
    $cplBase =~ tr/ACGTRYKMSWBDHV/TGCAYRMKWSVHDB/;
    $monoSubTotals{ substr( $sBase, 0, 1 ) } = 1
        if ( $monoSubTotals{ substr( $sBase, 0, 1 ) } == 0 );
    $monoSubTotals{ substr( $cplBase, 0, 1 ) } = 1
        if ( $monoSubTotals{ substr( $cplBase, 0, 1 ) } == 0 );
    $monoSubProb{ $sBase . " "
          . $cplBase } = (
       (
         ( $monoSubCounts{$sBase} / $monoSubTotals{ substr( $sBase, 0, 1 ) } ) +
             (
           $monoSubCounts{$cplBase} / $monoSubTotals{ substr( $cplBase, 0, 1 ) }
             )
       ) / 2
          );
  }

  return ( \%monoSubProb );
}
1;
