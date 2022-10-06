#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SeedAlignment.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An for holding a Dfam seed alignment dataset for a single
##      family.
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
# bless( {
#          'property' => value,
#           ...
#        }
#     }, 'SeedAlignment' );
#
###############################################################################
# ChangeLog
#
#     $Log: SeedAlignment.pm,v $
#     Revision 1.2  2017/04/05 00:03:31  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#

=head1 NAME

SeedAlignment

=head1 SYNOPSIS

use SeedAlignment

Usage: 

    my $familySeed = SeedAlignment->new();

=head1 DESCRIPTION

A class for storing the contents of a Dfam family seed alignment.

=head1 SEE ALSO

=over 4

SeedAlignmentCollection

=back

=head1 COPYRIGHT

Copyright 2017 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package SeedAlignment;
use strict;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use constant ExportComplete => 1;
use constant ExportHeadersOnly => 2;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "SeedAlignment";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new
{
  my $class          = shift;
  my %nameValuePairs = @_;
  my $subroutine     = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  # Allow import of values
  if ( %nameValuePairs )
  {
    while ( my ( $name, $value ) = each( %nameValuePairs ) )
    {

      # RMH: Perl optimisation, the calls to _ucFirst were
      #      costly.
      my $method = "set" . uc( substr( $name, 0, 1 ) ) . substr( $name, 1 );

      # RMH: Perl optimisation, the calls to ->can are really expensive
      #unless ( $this->can( $method ) ) {
      #  croak(
      #     "SeedAlignment::add: Instance variable $name doesn't exist." . "" );
      #}
      eval { $this->$method( $value ); };
      if ( $@ )
      {
        croak( "$subroutine: Instance variable $name doesn't exist." . "" );
      }
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##

=head2 clone()

  Use: my $newObj = $obj->clone();

  Clone a SeedAlignment *duplicating* all the values of the old
  object in the new one.

=cut

##-------------------------------------------------------------------------##
sub clone
{
  my $this = shift;

  my %newHash = %{$this};
  my $newObj  = \%newHash;

  bless $newObj, ref( $this );

  return $newObj;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setId()

  Use: my $value    = getId();
  Use: my $oldValue = setId( $value );

  Get/Set the id.

=cut

##-------------------------------------------------------------------------##
sub getId
{
  my $obj = shift;

  my $value = $obj->{'id'};

  return $value;
}

sub setId
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'id'};
  $obj->{'id'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setComment()

  Use: my $value    = getComment();
  Use: my $oldValue = setComment( $value );

  Get/Set the Comment.

=cut

##-------------------------------------------------------------------------##
sub getComment
{
  my $obj = shift;

  my $value = $obj->{'comment'};

  return $value;
}

sub setComment
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'comment'};
  $obj->{'comment'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setDescription()

  Use: my $value    = getDescription();
  Use: my $oldValue = setDescription( $value );

  Get/Set the description.

=cut

##-------------------------------------------------------------------------##
sub getDescription
{
  my $obj = shift;

  my $value = $obj->{'description'};

  return $value;
}

sub setDescription
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'description'};
  $obj->{'description'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setComments()

  Use: my $value    = getComments();
  Use: my $oldValue = setComments( $value );

  Get/Set the comments.

=cut

##-------------------------------------------------------------------------##
sub getComments
{
  my $obj = shift;

  my $value = $obj->{'comments'};

  return $value;
}

sub setComments
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'comments'};
  $obj->{'comments'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setCuratorComments()

  Use: my $value    = getCuratorComments();
  Use: my $oldValue = setCuratorComments( $value );

  Get/Set the curatorComments.

=cut

##-------------------------------------------------------------------------##
sub getCuratorComments
{
  my $obj = shift;

  my $value = $obj->{'curatorComments'};

  return $value;
}

sub setCuratorComments
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'curatorComments'};
  $obj->{'curatorComments'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setAlignmentMethod()

  Use: my $value    = getAlignmentMethod();
  Use: my $oldValue = setAlignmentMethod( $value );

  Get/Set the alignmentMethod.

=cut

##-------------------------------------------------------------------------##
sub getAlignmentMethod
{
  my $obj = shift;

  my $value = $obj->{'alignmentMethod'};

  return $value;
}

sub setAlignmentMethod
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'alignmentMethod'};
  $obj->{'alignmentMethod'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setClassification()

  Use: my $value    = getClassification();
  Use: my $oldValue = setClassification( $value );

  Get/Set the classification.

=cut

##-------------------------------------------------------------------------##
sub getClassification
{
  my $obj = shift;

  my $value = $obj->{'classification'};

  return $value;
}

sub setClassification
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'classification'};
  $obj->{'classification'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSeqCount()

  Use: my $value    = getSeqCount();
  Use: my $oldValue = setSeqCount( $value );

  Get/Set the seqCount.

=cut

##-------------------------------------------------------------------------##
sub getSeqCount
{
  my $obj = shift;

  my $value = $obj->{'seqCount'};

  return $value;
}

sub setSeqCount
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'seqCount'};
  $obj->{'seqCount'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setRfLine()

  Use: my $value    = getRfLine();
  Use: my $oldValue = setRfLine( $value );

  Get/Set the rfLine.

=cut

##-------------------------------------------------------------------------##
sub getRfLine
{
  my $obj = shift;

  my $value = $obj->{'rfLine'};

  return $value;
}

sub setRfLine
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'rfLine'};
  $obj->{'rfLine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setAuthors()

  Use: my $value    = getAuthors();
  Use: my $oldValue = setAuthors( $value );

  Get/Set the Authors.

=cut

##-------------------------------------------------------------------------##
sub getAuthors
{
  my $obj = shift;

  my $value = $obj->{'authors'};

  return $value;
}

sub setAuthors
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'authors'};
  $obj->{'authors'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setTSD()

  Use: my $value    = getTSD();
  Use: my $oldValue = setTSD( $value );

  Get/Set the TSD.

=cut

##-------------------------------------------------------------------------##
sub getTSD
{
  my $obj = shift;

  my $value = $obj->{'tsd'};

  return $value;
}

sub setTSD
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'tsd'};
  $obj->{'tsd'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setDatabaseRef()

  Use: my $value    = getDatabaseRef();
  Use: my $oldValue = setDatabaseRef( $value );

  Get/Set the DatabaseRef.

=cut

##-------------------------------------------------------------------------##
sub getDatabaseRef
{
  my $obj = shift;

  my $value = $obj->{'databaseRef'};

  return $value;
}

sub setDatabaseRef
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'databaseRef'};
  $obj->{'databaseRef'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setBuildMethod()

  Use: my $value    = getBuildMethod();
  Use: my $oldValue = setBuildMethod( $value );

  Get/Set the BuildMethod.

=cut

##-------------------------------------------------------------------------##
sub getBuildMethod
{
  my $obj = shift;

  my $value = $obj->{'buildMethod'};

  return $value;
}

sub setBuildMethod
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'buildMethod'};
  $obj->{'buildMethod'} = $value;

  return $oldValue;
}



##-------------------------------------------------------------------------##

=head2 setCitation() 

  Use: my ( $oldPmid, $oldTitle, $oldAuthor, $oldJournal ) = 
         setCitation( $index, $pmid, $title, $author, $journal );

  Set the citation at a given position in the list.

=cut

##-------------------------------------------------------------------------##
sub setCitation
{
  my $this       = shift;
  my $index      = shift;
  my $pmid       = shift;
  my $title      = shift;
  my $author     = shift;
  my $journal    = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'citations'} } );

  my ( $oldPmid, $oldTitle, $oldAuthor, $oldJournal ) =
      $this->getCitation( $index );

  $this->{'citations'}->[ $index ]->{'pmid'}    = $pmid;
  $this->{'citations'}->[ $index ]->{'title'}   = $title;
  $this->{'citations'}->[ $index ]->{'author'}  = $author;
  $this->{'citations'}->[ $index ]->{'journal'} = $journal;

  return ( $oldPmid, $oldTitle, $oldAuthor, $oldJournal );
}

##-------------------------------------------------------------------------##

=head2 getCitation() 

  Use: my ( $pmid, $title, $author, $journal ) = getCitation( $index );

  Get the citation at a given position in the list.

=cut

##-------------------------------------------------------------------------##
sub getCitation
{
  my $this       = shift;
  my $index      = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'citations'} } );

  return (
           $this->{'citations'}->[ $index ]->{'pmid'},
           $this->{'citations'}->[ $index ]->{'title'},
           $this->{'citations'}->[ $index ]->{'author'},
           $this->{'citations'}->[ $index ]->{'journal'}
  );
}

##-------------------------------------------------------------------------##

=head2 citationCount() 

  Use: my $count = citationCount();

  Return a count of the citation list.

=cut

##-------------------------------------------------------------------------##
sub citationCount
{
  my $this = shift;

  return $#{ $this->{'citations'} } + 1;
}

##-------------------------------------------------------------------------##

=head2 removeCitation() 

  Use: removeCitation( $index );

  Remove the citation at the given index.

=cut

##-------------------------------------------------------------------------##
sub removeCitation
{
  my $this       = shift;
  my $index      = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'citations'} } );

  splice( @{ $this->{'citations'} }, $index, 1 );
}

##-------------------------------------------------------------------------##

=head2 addCitation() 

  Use: addCitation( $pmid, $title, $author, $journal );

  Add a citation to the end of the list.

=cut

##-------------------------------------------------------------------------##
sub addCitation
{
  my $this    = shift;
  my $pmid    = shift;
  my $title   = shift;
  my $author  = shift;
  my $journal = shift;

  push @{ $this->{'citations'} },
      {
        pmid    => $pmid,
        title   => $title,
        author  => $author,
        journal => $journal
      };
}

##-------------------------------------------------------------------------##

=head2 clearCitations() 

  Use: clearCitations();

  Remove all citations from the list.

=cut

##-------------------------------------------------------------------------##
sub clearCitations
{
  my $this = shift;

  $this->{'citations'} = [];
}

##-------------------------------------------------------------------------##

=head2 setAlignment() 

  Use: my ( $oldAssemblyName, $oldSequenceName, $oldStart, 
            $oldEnd, $oldOrient, $oldSequence ) = 
         setAlignment( $index, $assemblyName, $sequenceName, 
                      $start, $end, $orient, $sequence );

  Set the alignment at a given position in the list.

=cut

##-------------------------------------------------------------------------##
sub setAlignment
{
  my $this         = shift;
  my $index        = shift;
  my $assemblyName = shift;
  my $sequenceName = shift;
  my $start        = shift;
  my $end          = shift;
  my $orient       = shift;
  my $sequence     = shift;
  my $subroutine   = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'alignments'} } );

  my @oldValues = $this->getAlignment( $index );

  $this->{'alignments'}->[ $index ]->{'assemblyName'} = $assemblyName;
  $this->{'alignments'}->[ $index ]->{'sequenceName'} = $sequenceName;
  $this->{'alignments'}->[ $index ]->{'start'}        = $start;
  $this->{'alignments'}->[ $index ]->{'end'}          = $end;
  $this->{'alignments'}->[ $index ]->{'orient'}       = $orient;
  $this->{'alignments'}->[ $index ]->{'sequence'}     = $sequence;

  return ( @oldValues );
}

##-------------------------------------------------------------------------##

=head2 getAlignment() 

  Use: my ( $assemblyName, $sequenceName, $start, $end, 
            $orient, $sequence ) = getAlignment( $index );

  Get the alignment at a given position in the list.

=cut

##-------------------------------------------------------------------------##
sub getAlignment
{
  my $this       = shift;
  my $index      = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'alignments'} } );

  return (
           $this->{'alignments'}->[ $index ]->{'assemblyName'},
           $this->{'alignments'}->[ $index ]->{'sequenceName'},
           $this->{'alignments'}->[ $index ]->{'start'},
           $this->{'alignments'}->[ $index ]->{'end'},
           $this->{'alignments'}->[ $index ]->{'orient'},
           $this->{'alignments'}->[ $index ]->{'sequence'}
  );
}

##-------------------------------------------------------------------------##

=head2 alignmentCount() 

  Use: my $count = alignmentCount();

  Return a count of the alignment list.

=cut

##-------------------------------------------------------------------------##
sub alignmentCount
{
  my $this = shift;

  return $#{ $this->{'alignments'} } + 1;
}

##-------------------------------------------------------------------------##

=head2 removeAlignment() 

  Use: removeAlignment( $index );

  Remove the alignment at the given index.

=cut

##-------------------------------------------------------------------------##
sub removeAlignment
{
  my $this       = shift;
  my $index      = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'alignments'} } );

  splice( @{ $this->{'alignments'} }, $index, 1 );
}

##-------------------------------------------------------------------------##

=head2 addAlignment() 

  Use: addAlignment( $assemblyName, $sequenceName, 
                     $start, $end, $orient, $sequence );

  Add a alignment to the end of the list.

=cut

##-------------------------------------------------------------------------##
sub addAlignment
{
  my $this         = shift;
  my $assemblyName = shift;
  my $sequenceName = shift;
  my $start        = shift;
  my $end          = shift;
  my $orient       = shift;
  my $sequence     = shift;

  push @{ $this->{'alignments'} },
      {
        assemblyName => $assemblyName,
        sequenceName => $sequenceName,
        start        => $start,
        end          => $end,
        orient       => $orient,
        sequence     => $sequence
      };
}

##-------------------------------------------------------------------------##

=head2 clearAlignments() 

  Use: clearAlignments();

  Remove all alignments from the list.

=cut

##-------------------------------------------------------------------------##
sub clearAlignments
{
  my $this = shift;

  $this->{'alignments'} = [];
}

##-------------------------------------------------------------------------##

=head2 setClade() 

  Use: my $oldCladeName = setClade( $index, $cladeName );

  Set the clade at a given position in the list.

=cut

##-------------------------------------------------------------------------##
sub setClade
{
  my $this       = shift;
  my $index      = shift;
  my $cladeName  = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'clades'} } );

  my $oldValue = $this->getClade( $index );

  $this->{'clades'}->[ $index ]->{'cladeName'} = $cladeName;

  return ( $oldValue );
}

##-------------------------------------------------------------------------##

=head2 getClade() 

  Use: my $cladeName = getClade( $index );

  Get the clade at a given position in the list.

=cut

##-------------------------------------------------------------------------##
sub getClade
{
  my $this       = shift;
  my $index      = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'clades'} } );

  return ( $this->{'clades'}->[ $index ]->{'cladeName'} );
}

##-------------------------------------------------------------------------##

=head2 cladeCount() 

  Use: my $count = cladeCount();

  Return a count of the clade list.

=cut

##-------------------------------------------------------------------------##
sub cladeCount
{
  my $this = shift;

  return $#{ $this->{'clades'} } + 1;
}

##-------------------------------------------------------------------------##

=head2 removeClade() 

  Use: removeClade( $index );

  Remove the clade at the given index.

=cut

##-------------------------------------------------------------------------##
sub removeClade
{
  my $this       = shift;
  my $index      = shift;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  croak "$subroutine: index $index is out of bounds."
      if ( $index < 0 || $index > $#{ $this->{'clades'} } );

  splice( @{ $this->{'clades'} }, $index, 1 );
}

##-------------------------------------------------------------------------##

=head2 addClade() 

  Use: addClade( $cladeName );

  Add a clade to the end of the list.

=cut

##-------------------------------------------------------------------------##
sub addClade
{
  my $this      = shift;
  my $cladeName = shift;

  push @{ $this->{'clades'} }, { cladeName => $cladeName };
}

##-------------------------------------------------------------------------##

=head2 clearClades() 

  Use: clearClades();

  Remove all clades from the list.

=cut

##-------------------------------------------------------------------------##
sub clearClades
{
  my $this = shift;

  $this->{'clades'} = [];
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2

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
     TD   TSD:                        Target Site Duplication String.
     DR   Database Reference:         Reference to external database. 
     CC   Long Description:           Multi-line description of the family.
     **   Curators Comments:          Dfam internal-use field for curator comments.

  #=GC
  
  Optional fields:
  ----------------
     RF        Reference annotation   Often the consensus DNA is used as a reference
                                      or simple "x" for match columns and "." for 
                                      insert columns.
  
=cut

##-------------------------------------------------------------------------##
sub read_stockholm
{
  my $this = shift;
  my $in   = shift;
  my $input;
  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  if ( ref $in eq "GLOB" )
  {
    $input = [ <$in> ];
  } else
  {
    $input = $in;
  }

  my $inCitation  = 0;
  my $pmid        = "";
  my $title       = "";
  my $author      = "";
  my $journal     = "";
  my $sequenceLen = 0;
  foreach ( @{$input} )
  {

    # Start of the data
    next if ( /^\# STOCKHOLM/ );

    # End of the data
    last if ( /^\/\// );

    # Specific GC line we care about ( RF )
    if ( /^\#=GC\s+RF\s+(\S+)/ )
    {
      $this->setRfLine( $1 );
      next;
    }

    # Save the first identifier line
    if ( /^\#=GF\s+ID\s+(\S+)/ && !$this->getId() )
    {

      # Only keep the first word in the field.
      $this->setId( $1 );
      next;
    }

    # Save the description lines
    if ( /^\#=GF\s+DE\s+(\S.*)$/ )
    {
      if ( $this->getDescription() )
      {
        $this->setDescription( $this->getDescription() . $1 . "\n" );
      } else
      {
        $this->setDescription( $1 . "\n" );
      }
      next;
    }
    
    # Save the comment lines
    if ( /^\#=GF\s+CC\s+(\S.*)$/ )
    {
      if ( $this->getComments() )
      {
        $this->setComments( $this->getComments() . $1 . "\n" );
      } else
      {
        $this->setComments( $1 . "\n" );
      }
      next;
    }

    # Save the curator comment lines
    if ( /^\#=GF\s+\*\*\s+(\S.*)$/ )
    {
      if ( $this->getCuratorComments() )
      {
        $this->setCuratorComments( $this->getCuratorComments() . $1 . "\n" );
      } else
      {
        $this->setCuratorComments( $1 . "\n" );
      }
      next;
    }


    # Save the classification of the entry
    if ( /^\#=GF\s+TP\s+(\S.*)$/ )
    {
      $this->setClassification( $1 );
      next;
    }
    
    # Save the DR entry
    if ( /^\#=GF\s+DR\s+(\S.*)$/ )
    {
      $this->setDatabaseRef( $1 );
      next;
    }
    
    # Save the AU entry
    if ( /^\#=GF\s+AU\s+(\S.*)$/ )
    {
      $this->setAuthors( $1 );
      next;
    }
     
    # Save the TD entry
    if ( /^\#=GF\s+TD\s+(\S.*)$/ )
    {
      $this->setTSD( $1 );
      next;
    }
    
    # Save the BM entry
    if ( /^\#=GF\s+BM\s+(\S.*)$/ )
    {
      $this->setBuildMethod( $1 );
      next;
    }

    # Save the clades
    if ( /^\#=GF\s+OC\s+(\S.*)$/ )
    {
      $this->addClade( $1 );
      next;
    }

    # Citations
    if ( /^\#=GF\s+RN\s+.*$/ )
    {
      if ( $inCitation )
      {

        # Save
        $this->addCitation( $pmid, $title, $author, $journal );
      }
      $pmid       = "";
      $title      = "";
      $author     = "";
      $journal    = "";
      $inCitation = 1;
      next;
    }
    if ( $inCitation && /^\#=GF\s+RM\s+(\S+)/ )
    {
      $pmid = $1;
      next;
    }
    if ( $inCitation && /^\#=GF\s+RT\s+(\S.*)$/ )
    {
      $title = $1;
      next;
    }
    if ( $inCitation && /^\#=GF\s+RA\s+(\S.*)$/ )
    {
      $author = $1;
      next;
    }
    if ( $inCitation && /^\#=GF\s+RL\s+(\S.*)$/ )
    {
      $journal = $1;
      next;
    }

    # Save unused GF lines for later parsing
    if ( /^\#=GF/ )
    {
      push( @{ $this->{'GF_lines'} }, $_ );
      next;
    }

    # Save sequence lines
    if ( /^([^\#]\S+)\s+([A-Za-z\.\-]+)\s*$/ )
    {
      my $origName = $1;
      my $sequence = $2;

      # Vaildate that all sequence lines are the same length
      my $tmpSeqLen = length( $sequence );
      if ( $sequenceLen > 0 && $tmpSeqLen != $sequenceLen )
      {
        croak "\n$subroutine:\n"
            . "  Line [$_] contains a sequence line that is\n"
            . "  not the same length as previous sequences.\n";

      }
      $sequenceLen = $tmpSeqLen;

      # Validate that all sequences only contain DNA or IUB codes.
      # Being strict here and not accepting 'x'.
      if ( $sequence !~ /^[AGCTYRWSKMDVHBNagctyrwskmdvhbn\.\-]+$/ )
      {
        croak "\n$subroutine:\n"
            . "  Line [$_] contains one or more non-DNA/IUB characters.\n";

      }

      my $assemblyName = "";
      my $sequenceName = "";
      my $start        = "";
      my $end          = "";
      my $orient       = "+";

      # Validate name
      if ( $origName =~ /^(\S+)\:(\S+)\:(\d+)-(\d+)(_[\+\-])?$/ )
      {
        $assemblyName = $1;
        $sequenceName = $2;
        $start        = $3;
        $end          = $4;
        my $optOrient = $5;
        if ( $optOrient ne "" ) {
          if ( $optOrient =~ /_([\+\-])/ ) {
            $orient = $1;
          }else {
            croak "\n$subroutine:\n"
                . "  Line [$origName ...] contains an unparable identifier!\n";
          }
        }elsif ( $end < $start )
        {
          $end    = $3;
          $start  = $4;
          $orient = "-";
        }
      } elsif ( $origName =~ /^(\S+)\:(\d+)-(\d+)(_[\+\-])?$/ )
      {
        $sequenceName = $1;
        $start        = $2;
        $end          = $3;
        my $optOrient = $4;
        if ( $optOrient ne "" ) {
          if ( $optOrient =~ /_([\+\-])/ ) {
            $orient = $1;
          }else {
            croak "\n$subroutine:\n"
                . "  Line [$origName ...] contains an unparable identifier!\n";
          }
        }elsif ( $end < $start )
        {
          $end    = $2;
          $start  = $3;
          $orient = "-";
        }
      } else
      {
        # This message needs to be extremely helpful so we encourage users
        # to user standard naming conventions.
        croak "\n$subroutine:\n"
            . "  Line [$origName ...] contains a sequence\n"
            . "  name that does not conform to the Dfam standard.  Please name the\n"
            . "  sequences using either the \"assemblyID:sequenceID:start-end\", \n"
            . "  \"assemblyID:sequenceID:start-end_orient\", or \n"
            . "  \"sequenceID:start-end\" format ( 1-based coordinates ).  Also\n"
            . "  note that sequences on the reverse strand are denoted by\n"
            . "  either including the \"_-\" suffix or by reversing the start and\n"
            . "  end coordinates (ie. \"hg38:chr1:1103-500\").\n"
            . "  Please use public identifiers for assembly/sequence wherever\n"
            . "  possible.\n\n";
      }

      $this->addAlignment( $assemblyName, $sequenceName, $start, $end, $orient,
                           $sequence );

      next;
    }

    # Blank line? fine. Comment? fine. Anything else? Forget it
    croak "$subroutine: Line [$_] is not valid stockholm format"
        if ( /\w/ && !/^\#/ );
  }
  if ( $inCitation )
  {
    $this->addCitation( $pmid, $title, $author, $journal );
  }

}


##-------------------------------------------------------------------------##

=head2

  Use: toString( $format );

  Create stockholm string representation of the object.

      $format        : SeedAlignment::ExportComplete      [ default ]
                       SeedAlignment::ExportHeadersOnly 

=cut

##-------------------------------------------------------------------------##
sub toString
{
  my $obj           = shift;
  my $format        = shift;
  my $displayParams = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  $format = SeedAlignment::ExportComplete
      if ( !defined $format );

  my $retStr = "";

  $retStr .= "# STOCKHOLM 1.0\n";
  $retStr .= "#=GF ID    " . $obj->getId() . "\n";
  if ( $obj->getDescription() )
  {
    my @lines = split( /[\n\r]+/, $obj->getDescription() );
    foreach my $line ( @lines )
    {
      $retStr .= "#=GF DE    $line\n";
    }
  }
  if ( $obj->getAuthors() ) {
    $retStr .= "#=GF AU    " . $obj->getAuthors() . "\n";
  }
  if ( $obj->getComment() )
  {
    my @lines = split( /[\n\r]+/, $obj->getComment() );
    foreach my $line ( @lines )
    {
      $retStr .= "#=GF CC    $line\n";
    }
  }
  $retStr .= "#=GF TP    " . $obj->getClassification() . "\n"
      if ( $obj->getClassification() );
  for ( my $i = 0 ; $i < $obj->cladeCount() ; $i++ )
  {
    $retStr .= "#=GF OC    " . $obj->getClade( $i ) . "\n";
  }
  if ( $obj->getTSD() ) {
    $retStr .= "#=GF TD    " . $obj->getTSD() . "\n";
  }
  for ( my $i = 0 ; $i < $obj->citationCount() ; $i++ )
  {
    my ( $pmid, $title, $author, $journal ) = $obj->getCitation( $i );
    $retStr .= "#=GF RN    [" . ( $i + 1 ) . "]\n";
    $retStr .= "#=GF RM    $pmid\n";
    $retStr .= "#=GF RT    $title\n" if ( $title );
    $retStr .= "#=GF RA    $author\n" if ( $author );
    $retStr .= "#=GF RL    $journal\n" if ( $journal );
  }
  if ( $obj->getDatabaseRef() ) {
    $retStr .= "#=GF DR    " . $obj->getDatabaseRef() . "\n";
  }
  if ( $obj->getComments() )
  {
    my @lines = split( /[\n\r]+/, $obj->getComments() );
    foreach my $line ( @lines )
    {
      $retStr .= "#=GF CC    $line\n";
    }
  }
  if ( $obj->getCuratorComments() )
  {
    my @lines = split( /[\n\r]+/, $obj->getCuratorComments() );
    foreach my $line ( @lines )
    {
      $retStr .= "#=GF **    $line\n";
    }
  }
  if ( $obj->getBuildMethod() ) {
    $retStr .= "#=GF BM    " . $obj->getBuildMethod() . "\n";
  }
  if ( $format != SeedAlignment::ExportHeadersOnly ) { 
    # These two "headers" are special as they relate directly to the MSA data.  
    # So..let's consider these as part of the sequence data.
    $retStr .= "#=GF SQ    " . $obj->alignmentCount() . "\n";
    $retStr .= "#=GC RF    " . $obj->getRfLine() . "\n"
        if ( $obj->getRfLine() );
  
    for ( my $i = 0 ; $i < $obj->alignmentCount() ; $i++ )
    {
      my ( $assemblyName, $sequenceName, $start, $end, $orient, $sequence ) =
          $obj->getAlignment( $i );
  
      my $id;
      $id = "$assemblyName:" if ( $assemblyName );
      $id .= "$sequenceName:";
      if ( $orient eq "+" )
      {
        $id .= "$start-$end";
      } else
      {
        $id .= "$end-$start";
      }
      $retStr .= "$id    $sequence\n";
    }
    $retStr .= "//\n";
  }
  return ( $retStr );
}

sub reverseComplementAlignment
{
  my $obj = shift;

  my $rfLine = $obj->getRfLine();
  if ( $rfLine ne "" )
  {
    $rfLine = reverse( $rfLine );

    # In case a consensus sequence is used in the RF line
    # make sure we complement it.
    $rfLine =~ tr/ACGTRYWSKMNXBDHV/TGCAYRSWMKNXVHDB/;
    $rfLine =~ tr/acgtrywskmnxbdhv/tgcayrswmknxvhdb/;
    $obj->setRfLine( $rfLine );
  }

  for ( my $i = 0 ; $i < $obj->alignmentCount() ; $i++ )
  {
    my ( $assemblyName, $sequenceName, $start, $end, $orient, $sequence ) =
        $obj->getAlignment( $i );

    if ( $orient eq "+" )
    {
      $orient = "-";
    } else
    {
      $orient = "+";
    }

    # Reverse the sequence
    $sequence = reverse( uc( $sequence ) );
    $sequence =~ tr/ACGTRYWSKMNXBDHV/TGCAYRSWMKNXVHDB/;

    $obj->setAlignment( $i, $assemblyName, $sequenceName, $start, $end, $orient,
                        $sequence );
  }

}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = _toString([$this]);
##
##      $this         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
sub _toString
{
  my $this = shift;
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

1;
