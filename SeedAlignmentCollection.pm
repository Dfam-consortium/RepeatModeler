#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SeedAlignmentCollection.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An datastructure for holding a set of SeedAlignment
##      objects.  This object can read/write Stockholm files
##      in the Dfam dialect.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2017 Developed by
#* Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless( { 'seed_objects' => [ $SeedObjectRef1,
#                              $SeedObjectRef2,
#                              ..
#                            ]
#        }
#     }, 'SeedAlignmentCollection' );
#
###############################################################################
# ChangeLog
#
#     $Log: SeedAlignmentCollection.pm,v $
#     Revision 1.2  2017/04/05 00:03:31  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#
#

=head1 NAME

SeedAlignmentCollection

=head1 SYNOPSIS

use SeedAlignmentCollection

Usage: 

    my $stockholmFile = SeedAlignmentCollection->new();
    my $IN;
    open $IN,"<library.stk" or die;
    $stockholmFile->read_stockholm( $IN );

=head1 DESCRIPTION

An datastructure for holding a set of SeedAlignment objects.  This 
object can read/write Stockholm files in the Dfam dialect.

=head1 INSTANCE METHODS

=cut 

package SeedAlignmentCollection;
use strict;
use SeedAlignment;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
require ArrayList;
require Exporter;

@ISA = qw(ArrayList Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "SeedAlignmentCollection";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  my $this = $class->SUPER::new( @_ );

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 read_stockholm()

  Use: $obj->read_stockholm(\*INPUT);
    or
       $obj->read_stockholm(\@array_ref);

  Reads the Dfam dialect of the Stockholm format.

  Minimal Example
  ===============

  # STOCKHOLM 1.0
  #=GF ID    Jumbo1
  #=GF DE    A very common DNA transposon with two deletion products Jumbo1A
  #=GF DE    and Jumbo1B. 
  #=GF SE    Predicted; RepeatModeler
  #=GF TP    Interspersed_Repeat;Unknown
  #=GF OC    Colius striatus
  #=GF SQ    2
  #=GC RF                         xxxx.xxx.x.x...x.x.x
  colStr1:KK551537:690782-691123  .GCT.GGGAT.G...CTGCG
  colStr1:KK551537:1038-2834      AGCTTGGG.TTGTACC.G.G
  //

  Complex Example
  ===============

  # STOCKHOLM 1.0
  #=GF ID    Jumbo1
  #=GF DE    A very common DNA transposon with two deletion products Jumbo1A
  #=GF DE    and Jumbo1B. 
  #=GF SE    Predicted; RepeatModeler
  #=GF TP    Interspersed_Repeat;Transposable_Element;DNA;TIR;hAT;Tip100
  #=GF OC    Eutheria
  #=GF OC    Metatheria
  #=GF RN    [1]
  #=GF RM    382832
  #=GF RN    [2]
  #=GF RM    28834
  #=GF DR    DPTEdb; af31b_dna; 
  #=GF SQ    2
  #=GC RF                         xxxx.xxx.x.x...x.x.x
  hg38:chr12:690782-691123        .GCT.GGGAT.G...CTGCG
  hg38:chr2:38282-38399           AGCTTGGG.TTGTACC.G.G
  //



  Recommended features
  ====================

  #=GF

  Compulsory fields:
  ------------------
     ID   Identification:             One word name for family.
     DE   Definition:                 Short description of family.
     AU   Author:                     Authors of the entry.
     SE   Source of seed:             The source suggesting the seed members 
                                      belong to one family.
     TP   Type/Classification:        Classifcation of family -- Presently we use the
                                      new unified RepeatMasker classification 
                                      heirarchy.
     OC   Clade:                      Organism (clade, etc.) Multiple OC records are
                                      allowed.
     SQ   Sequence:                   Number of sequences in alignment.


  Optional fields:
  ----------------
     RN   Reference Number:           Reference Number.
     RM   Reference PubMed:           Pubmed reference number.
     RT   Reference Title:            Reference Title. 
     RA   Reference Author:           Reference Author
     RL   Reference Location:         Journal location. 
     DR   Database Reference:         Reference to external database. 


  #=GC

  Optional fields:
  ----------------
     RF        Reference annotation   Often the consensus DNA is used as a reference
                                      or simple "x" for match columns and "." for 
                                      insert columns.

=cut

##-------------------------------------------------------------------------##
sub read_stockholm {
  my $this = shift;
  my $in   = shift;
  my $input;

  if ( ref $in eq "GLOB" ) {
    $input = [ <$in> ];
  }
  else {
    $input = $in;
  }

  my @data;
  my $inSection = 0;
  foreach ( @$input ) {
    if ( /^\#\s+STOCKHOLM/ ) {
      if ( $inSection ) {
        croak $CLASS
            . "::read_stockholm: Two headings found without an intervening '//' line.";
      }
      $inSection = 1;
    }

    if ( /^\/\// && $inSection && @data ) {

      # Process it
      push @data, $_;
      my $seedAlignObj = SeedAlignment->new();
      $seedAlignObj->read_stockholm( \@data );
      $this->add( $seedAlignObj );
      @data      = ();
      $inSection = 0;
    }

    push @data, $_
        if ( $inSection );
  }
  if ( @data ) {
    croak $CLASS
        . "::read_stockholm: Incomplete Stockholm section at the end of the file.";
  }
}

##-------------------------------------------------------------------------##

=head2 write()

  Use: $obj->write( $file );

    $file : String

=cut

##-------------------------------------------------------------------------##
sub write {
  my $this     = shift;
  my $filename = shift;

  open OUT, ">$filename"
      or croak "SeedAlignmentCollection::write - "
      . "Error...could not open file $filename "
      . "for writing!\n";

  for ( my $i = 0 ; $i < $this->size() ; $i++ ) {
    print OUT $this->get( $i )->toString();
  }
  close OUT;
}

##-------------------------------------------------------------------------##

=head2 toString()

  Use: my $str = $obj->toString();

  Create a monolithic string with all the object data
  formatted for output.

=cut

##-------------------------------------------------------------------------##
sub toString() {
  my $this          = shift;
  my $alignmentMode = shift;

  my $str = "";
  for ( my $i = 0 ; $i < $this->size() ; $i++ ) {
    $str .= $this->get( $i )->toString();
  }
  return $str;
}

1;
