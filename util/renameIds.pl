#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) renameIDs
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A simple script to assist in renaming family identifiers
##      in a Stockholm file.
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
#     $Log: renameIds.pl,v $
#     Revision 1.2  2017/04/05 00:03:32  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#      - Create a translation log for old names
#      - Allow custom format syntax
#

=head1 NAME

renameIDs - Rename identifiers in a Stockholm file

=head1 SYNOPSIS

  renameIDs.pl [-version] [-format <format_string>] 
               [-trans <translation_file>]
               [-old <old_id> -new <new_id>] 
               <input_stockholm_file> <output_stockholm_file>

  Rename a single identifier in a stockholm file:

    renameIDs.pl -old <old_id> -new <new_id> <input_stockholm_file>
                 <output_stockholm_file>

    ie.
        renameIDs.pl -old rnd-1_family-2 -new Piper-1a myfile.stk 
                     outfile.stk

  Rename all identifiers in the file using the given format string.
  A format string may contain any combination of alpha numeric
  characters ( [a-zA-Z0-9] ), symbols "-", "_", or tags.
  Tags are two character specifiers for data elements which can
  be inserted in the identifier by the script.  The following
  tags are currently supported:

      %i  - Index starting at 1 and incrementing for 
            each duplicate identifier.

      %t  - RepeatMasker type for the assigned classification
            or "Unk" if not defined.

      %s  - RepeatMasker subtype for assigned classification
            or "Unk" if not defined.

  For example, to create id's like: 
    <class_abbreviation>-<index_prefix>#_hum
    
  The format specifier would be: 
  
        renameIDs.pl -format "%t-%i_hum" myfile.stk outfile.stk

  Rename all identifiers using a translation table provided in
  the form of a tab separated text file.  Each line should be
  of the form: <old_name><tab><new_name>

    renameIDs.pl -trans <translation.tsv> <input_stockholm_file>
                 <output_stockholm_file>
  

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2017 Robert Hubley, Institute for Systems Biology

=head1 LICENSE

This code may be used in accordance with the Open Source License v. 3.0
http://opensource.org/licenses/OSL-3.0

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
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/..";
use SeedAlignmentCollection;
use SeedAlignment;

#
# Version
#  -- NOTE: This is filled in by configure
my $Version = "open-1.0.11";
$Version = "DEV" if ( $Version =~ /\#VERSION\#/ );

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
                    '-version',    # print out the version and exit
                    '-old=s',
                    '-new=s',
                    '-trans=s',
                    '-format=s',
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

#
# ARGV Processing
#
if ( !@ARGV )
{
  usage();
}

my $sFilename     = $ARGV[ 0 ];
my $outFile       = $ARGV[ 1 ];
my $stockholmFile = SeedAlignmentCollection->new();
open my $IN, "<$sFilename"
    or die "Could not open up stockholm file $sFilename for reading!\n";
$stockholmFile->read_stockholm( $IN );
close $IN;

my %translationTable;
if ( $options{'trans'} )
{
  open IN, "<$options{'trans'}"
      or die "Could not open translation table: $options{'trans'}\n";
  while ( <IN> )
  {
    if ( /^\s*(\S+)\s*\t\s*(\S+)\s*$/ )
    {
      $translationTable{$1} = $2;
    } else
    {
      warn "Could not parse translation from: $_\n";
    }
  }
  close IN;
}

open OUT, ">$outFile" or die "Could not open $outFile for writing!\n";

my $idxPrefix = "";
$idxPrefix = $options{'idxpre'} if ( $options{'idxpre'} );
my %idxs   = ();
my %unique = ();
for ( my $i = 0 ; $i < $stockholmFile->size() ; $i++ )
{
  my $seedAlign = $stockholmFile->get( $i );
  my $id        = $seedAlign->getId();
  if ( $options{'old'} && $id =~ /$options{'old'}/i )
  {
    $seedAlign->setId( $options{'new'} );
  } elsif ( $options{'trans'} )
  {
    if ( defined $translationTable{$id} )
    {
      $seedAlign->setId( $translationTable{$id} );
    }
  } elsif ( $options{'format'} )
  {
    my $format = $options{'format'};

    # Sanity check format
    if ( $format !~ /^[a-zA-Z0-9\-\_\%]+$/ )
    {
      print "\nFormat \"$format\" contains invalid characters.\n\n";
      exit;
    }

    my $class = $seedAlign->getClassification();

    my $type    = "Unk";
    my $subtype = "Unk";
    if ( defined $class && $class ne "" )
    {
      ( $type, $subtype ) = getTypeSubtype( $class );
    }

    my $newID = $format;
    $newID =~ s/\%t/$type/g;
    $newID =~ s/\%s/$subtype/g;

    if ( $newID =~ /\%i/ )
    {
      if ( !exists $idxs{$newID} )
      {
        $idxs{$newID} = 0;
      }
      $idxs{$newID}++;
      $newID =~ s/\%i/$idxs{$newID}/g;
    }

    if ( $newID =~ /\%/ )
    {
      print "Error in format.  Unmatched \"\%\" specifier.\n";
      exit;
    }

    if ( exists $unique{$newID} )
    {
      print "Duplicate identifier \"$newID\"!  Try using a \"\%i\"\n"
          . "in your format string.\n";
      exit;
    }
    $unique{$newID} = 1;
    $seedAlign->setId( $newID );

    # TODO: Add alias
    # setDatabaseReference()...  need to work on SeedAlignment.pm
  }

  print OUT "" . $seedAlign->toString() . "\n";
}
close OUT;

sub getTypeSubtype
{
  my $class = shift;

  if ( $class !~ /^root;/ )
  {
    $class = "root;" . $class;
  }
  my $cache = {
    "root"                                           => [],
    "root;accidental"                                => [],
    "root;accidental;normally_non-integrating_virus" =>
        [ "Other", "DNA_virus" ],
    "root;other"                                    => [ "Other" ],
    "root;interspersed_repeat"                      => [],
    "root;interspersed_repeat;transposable_element" => [],
    "root;interspersed_repeat;transposable_element;dna_transposon" => [],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton"
        => [ "DNA", "Crypton" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-s"
        => [ "DNA", "Crypton-S" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-r"
        => [ "DNA", "Crypton-R" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-v"
        => [ "DNA", "Crypton-V" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-f"
        => [ "DNA", "Crypton-F" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-c"
        => [ "DNA", "Crypton-C" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-i"
        => [ "DNA", "Crypton-I" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-x"
        => [ "DNA", "Crypton-X" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-a"
        => [ "DNA", "Crypton-A" ],
"root;interspersed_repeat;transposable_element;dna_transposon;circular_dsdna_intermediate;crypton;crypton-h"
        => [ "DNA", "Crypton-H" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat"
        => [ "DNA" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;p_element"
        => [ "DNA", "P" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;p_element;fungi-specific_branch"
        => [ "DNA", "P-Fungi" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;kolobok_group"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;kolobok_group;hydra-specific_branch"
        => [ "DNA", "Kolobok-Hydra" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;kolobok_group;kolobok-h"
        => [ "DNA", "Kolobok-H" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;kolobok_group;kolobok-e"
        => [ "DNA", "Kolobok-E" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;kolobok_group;kolobok"
        => [ "DNA", "Kolobok" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;kolobok_group;t2"
        => [ "DNA", "Kolobok-T2" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element"
        => [ "DNA", "hAT" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;restless"
        => [ "DNA", "hAT-Restless" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;blackjack"
        => [ "DNA", "hAT-Blackjack" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hobo"
        => [ "DNA", "hAT-hobo" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;tag1"
        => [ "DNA", "hAT-Tag1" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hatw"
        => [ "DNA", "hAT-hATw" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;tip100"
        => [ "DNA", "hAT-Tip100" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hat6"
        => [ "DNA", "hAT-hAT6" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hat19"
        => [ "DNA", "hAT-hAT19" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hat1"
        => [ "DNA", "hAT-hAT1" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hatx"
        => [ "DNA", "hAT-hATx" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;pegasus"
        => [ "DNA", "hAT-Pegasus" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hatm"
        => [ "DNA", "hAT-hATm" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;hat5"
        => [ "DNA", "hAT-hAT5" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;activator"
        => [ "DNA", "hAT-Ac" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;hat_element;charlie"
        => [ "DNA", "hAT-Charlie" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;ginger"
        => [ "DNA", "Ginger" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element;transib"
        => [ "DNA", "CMC-Transib" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element;cmc"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element;cmc;enspm"
        => [ "DNA", "CMC-EnSpm" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element;cmc;chapaev_group"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element;cmc;chapaev_group;chapaev"
        => [ "DNA", "CMC-Chapaev" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element;cmc;chapaev_group;chapaev-3"
        => [ "DNA", "CMC-Chapaev-3" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;cacta_element;cmc;mirage"
        => [ "DNA", "CMC-Mirage" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;sola-group"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;sola-group;sola-1"
        => [ "DNA", "Sola-1" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;sola-group;sola-3"
        => [ "DNA", "Sola-3" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;sola-group;sola-2"
        => [ "DNA", "Sola-2" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;zisupton"
        => [ "DNA", "Zisupton" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;zator"
        => [ "DNA", "Zator" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;pif-like_elements"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;pif-like_elements;isl2eu"
        => [ "DNA", "PIF-ISL2EU" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;pif-like_elements;harbinger"
        => [ "DNA", "PIF-Harbinger" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;pif-like_elements;spy"
        => [ "DNA", "PIF-Spy" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;pif-like_elements;harbs"
        => [ "DNA", "PIF-HarbS" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;academ_group"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;academ_group;academ-h"
        => [ "DNA", "Academ-H" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;academ_group;academ-2"
        => [ "DNA", "Academ-2" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;academ_group;academ-1"
        => [ "DNA", "Academ-1" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;merlin"
        => [ "DNA", "Merlin" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;piggybac-like_element"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;piggybac-like_element;piggybac-a"
        => [ "DNA", "PiggyBac-A" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;piggybac-like_element;piggybac"
        => [ "DNA", "PiggyBac" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;piggybac-like_element;piggybac-x"
        => [ "DNA", "PiggyBac-X" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;mutator-like_element"
        => [ "DNA", "MULE" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;mutator-like_element;f"
        => [ "DNA", "MULE-F" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;mutator-like_element;mudr"
        => [ "DNA", "MULE-MuDR" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;mutator-like_element;nof"
        => [ "DNA", "MULE-NOF" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;dada"
        => [ "DNA", "Dada" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;novosib"
        => [ "DNA", "Novosib" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;dna_pol"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;dna_pol;casposon"
        => [ "DNA", "Casposons" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element"
        => [ "DNA", "TcMar" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;pogo-group"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;pogo-group;tigger"
        => [ "DNA", "TcMar-Tigger" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;pogo-group;fot1"
        => [ "DNA", "TcMar-Fot1" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;pogo-group;pogo"
        => [ "DNA", "TcMar-Pogo" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;tc2"
        => [ "DNA", "TcMar-Tc2" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;tc4"
        => [ "DNA", "TcMar-Tc4" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;mogwai"
        => [ "DNA", "TcMar-Mogwai" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;stowaway"
        => [ "DNA", "TcMar-Stowaway" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;mariner"
        => [ "DNA", "TcMar-Mariner" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;cweed"
        => [ "DNA", "TcMar-Cweed" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;sagan"
        => [ "DNA", "TcMar-Sagan" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;gizmo"
        => [ "DNA", "TcMar-Gizmo" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;ant1"
        => [ "DNA", "TcMar-Ant1" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;tc1"
        => [ "DNA", "TcMar-Tc1" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;isrm11"
        => [ "DNA", "TcMar-ISRm11" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;tc1-mariner-like_element;m44"
        => [ "DNA", "TcMar-m44" ],
"root;interspersed_repeat;transposable_element;dna_transposon;terminal_inverted_repeat;is3eu"
        => [ "DNA", "IS3EU" ],
"root;interspersed_repeat;transposable_element;dna_transposon;dna_polymerase"
        => [],
"root;interspersed_repeat;transposable_element;dna_transposon;dna_polymerase;maverick"
        => [ "DNA", "Maverick" ],
"root;interspersed_repeat;transposable_element;dna_transposon;rolling_circle"
        => [ "RC" ],
"root;interspersed_repeat;transposable_element;dna_transposon;rolling_circle;helitron-2"
        => [ "RC", "Helitron-2" ],
"root;interspersed_repeat;transposable_element;dna_transposon;rolling_circle;helitron-1"
        => [ "RC", "Helitron" ],
    "root;interspersed_repeat;transposable_element;retrotransposed_element" =>
        [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter"
        => [ "Retroposon" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter;r4-derived"
        => [ "SINE", "Dong-R4" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter;l1-dependent"
        => [ "Retroposon", "L1-dep" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter;l1-dependent;sva"
        => [ "Retroposon", "SVA" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter;l1-derived"
        => [ "SINE", "L1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter;l2-derived"
        => [ "SINE", "L2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter;rte-derived"
        => [ "Retroposon", "RTE-Derived" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;lacking_small_rna_pol_iii_promoter;i-derived"
        => [ "SINE", "I" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine"
        => [ "SINE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_5s_rna"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_5s_rna;no_or_unknown_core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_5s_rna;no_or_unknown_core;unknown_line-dependent"
        => [ "SINE", "tRNA-5S" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;7sl-rna_promoter"
        => [ "SINE", "7SL" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;7sl-rna_promoter;no-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;7sl-rna_promoter;no-core;l1-dependent"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;7sl-rna_promoter;no-core;l1-dependent;alu"
        => [ "SINE", "Alu" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;7sl-rna_promoter;no-core;l1-dependent;b2"
        => [ "SINE", "B2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;unknown_promoter"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;unknown_promoter;ceph-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;unknown_promoter;ceph-core;rte-end"
        => [ "SINE", "Ceph" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;unknown_promoter;mir-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;unknown_promoter;mir-core;rte-end"
        => [ "SINE", "Core-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;unknown_promoter;mir-core;unknown_line-dependent"
        => [ "SINE", "Core" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;v-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;v-core;cr1-end"
        => [ "SINE", "tRNA-V-CR1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;v-core;unknown_line-dependent"
        => [ "SINE", "tRNA-V" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;sauria-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;sauria-core;l2-end"
        => [ "SINE", "tRNA-Sauria-L2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;sauria-core;rte-end"
        => [ "SINE", "tRNA-Sauria-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;sauria-core;unknown_line-dependent"
        => [ "SINE", "tRNA-Sauria" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;ceph-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;ceph-core;rte-end"
        => [ "SINE", "tRNA-Ceph-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;ceph-core;unknown_line-dependent"
        => [ "SINE", "tRNA-Ceph" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;l1-dependent"
        => [ "SINE", "tRNA-L1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;i-end"
        => [ "SINE", "tRNA-I" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;l2-end"
        => [ "SINE", "tRNA-L2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;rex-end"
        => [ "SINE", "tRNA-Rex" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;cr1-end"
        => [ "SINE", "tRNA-CR1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;jockey-end"
        => [ "SINE", "tRNA-Jockey" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;r1-end"
        => [ "SINE", "R1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;tad1_end"
        => [ "SINE", "tRNA-Tad1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;rte-end"
        => [ "SINE", "tRNA-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;r2-end"
        => [ "SINE", "tRNA-R2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;unknown_line-dependent"
        => [ "SINE", "tRNA" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no_or_unknown_core;bovb-end"
        => [ "SINE", "RTE-BovB" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;deu-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;deu-core;i-end"
        => [ "SINE", "tRNA-Deu-I" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;deu-core;l2-end"
        => [ "SINE", "tRNA-Deu-L2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;deu-core;rte-end"
        => [ "SINE", "tRNA-Deu-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;deu-core;cr1-end"
        => [ "SINE", "tRNA-Deu-CR1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;deu-core;unknown_line-dependent"
        => [ "SINE", "tRNA-Deu" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;no-core;l1-dependent"
        => [ "SINE", "ID" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;v_and_mir-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;v_and_mir-core;l2-end"
        => [ "SINE", "tRNA-V-Core-L2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;meta-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;meta-core;unknown_line-dependent"
        => [ "SINE", "tRNA-Meta" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;mir-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;mir-core;l2-end"
        => [ "SINE", "MIR" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;mir-core;mermaid"
        => [ "SINE", "tRNA-Mermaid" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;mir-core;rte-end"
        => [ "SINE", "tRNA-Core-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_promoter;mir-core;unknown_line-dependent"
        => [ "SINE", "tRNA-Core" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter"
        => [ "SINE", "5S" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;sauria-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;sauria-core;rte-end"
        => [ "SINE", "5S-Sauria-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;no_or_unknown_core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;no_or_unknown_core;rte-end"
        => [ "SINE", "5S-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;deu-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;deu-core;l2-end"
        => [ "SINE", "5S-Deu-L2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;deu-core;unknown_line-dependent"
        => [ "SINE", "5S-Deu" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;mir-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;5s-rna_promoter;mir-core;rte-end"
        => [ "SINE", "5S-Core-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;u-rna_promoter"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;u-rna_promoter;no_or_unknown_core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;u-rna_promoter;no_or_unknown_core;unknown_line-dependent"
        => [ "SINE", "U" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_7sl_rna"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_7sl_rna;no-core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_7sl_rna;no-core;l1-dependent"
        => [ "SINE", "B4" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_7sl_rna;no_or_unknown_core"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line-dependent_retroposon;sine;trna_and_7sl_rna;no_or_unknown_core;unknown_line-dependent"
        => [ "SINE", "tRNA-7SL" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;penelope-like_elements"
        => [ "LINE", "Penelope" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element"
        => [ "LTR" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;trim"
        => [ "LTR", "TRIM" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;bel-pao"
        => [ "LTR", "Pao" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;ty1-copia"
        => [ "LTR", "Copia" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae"
        => [ "LTR", "Cassandra" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;spumaretrovirinae"
        => [ "LTR", "ERV-Foamy" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae"
        => [ "LTR", "ERV" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2;lenti"
        => [ "LTR", "ERV-Lenti" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2-group"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2-group;erv4"
        => [ "LTR", "ERV4" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2-group;erv2"
        => [ "LTR", "ERVK" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2-group;erv3"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2-group;erv3;ervl"
        => [ "LTR", "ERVL" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv2-group;erv3;malr"
        => [ "LTR", "ERVL-MaLR" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;retroviridae;orthoretrovirinae;erv1"
        => [ "LTR", "ERV1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;pararetroviridae"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;pararetroviridae;caulimoviridae"
        => [ "LTR", "Caulimovirus" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;long_terminal_repeat_element;gypsy-erv;gypsy"
        => [ "LTR", "Gypsy" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;tate"
        => [ "Unknown", "TATE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;inverted_long_terminal_repeat_elements"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;inverted_long_terminal_repeat_elements;tyrosine_recombinase_elements"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;inverted_long_terminal_repeat_elements;tyrosine_recombinase_elements;dirs"
        => [ "LTR", "DIRS" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;inverted_long_terminal_repeat_elements;tyrosine_recombinase_elements;viper"
        => [ "LTR", "Viper" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;retrotransposon;inverted_long_terminal_repeat_elements;tyrosine_recombinase_elements;ngaro"
        => [ "LTR", "Ngaro" ],
    "root;interspersed_repeat;transposable_element;retrotransposed_element;line"
        => [ "LINE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-i"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-i;odin"
        => [ "LINE", "CRE-Odin" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-i;cre"
        => [ "LINE", "CRE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-i;cre;cre-2"
        => [ "LINE", "CRE-2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-i;cre;cre-1"
        => [ "LINE", "CRE-1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-i;ambal"
        => [ "LINE", "CRE-Ambal" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;genie"
        => [ "LINE", "Genie" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r2-like"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r2-like;hero"
        => [ "LINE", "R2-Hero" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r2-like;nesl"
        => [ "LINE", "R2-NeSL" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r2-like;r2"
        => [ "LINE", "R2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;l1-like"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;l1-like;dre"
        => [ "LINE", "L1-DRE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;l1-like;zorro"
        => [ "LINE", "L1-Zorro" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;l1-like;l1-group"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;l1-like;l1-group;tx1"
        => [ "LINE", "L1-Tx1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;l1-like;l1-group;l1"
        => [ "LINE", "L1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;proto-1"
        => [ "LINE", "Proto1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r4-like"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r4-like;dualen"
        => [ "LINE", "Dualen" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r4-like;dong-r4"
        => [ "LINE", "Dong-R4" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-1;r4-like;deceiver"
        => [ "LINE", "Deceiver" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;proto-2"
        => [ "LINE", "Proto2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;rte-like"
        => [ "LINE", "RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;rte-like;orte"
        => [ "LINE", "RTE-ORTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;rte-like;rte-group"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;rte-like;rte-group;rte-x"
        => [ "LINE", "RTE-X" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;rte-like;rte-group;rte"
        => [ "LINE", "RTE-RTE" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;rte-like;rte-group;bovb"
        => [ "LINE", "RTE-BovB" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;tad1"
        => [ "LINE", "Tad1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;cr1-group"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;cr1-group;l2-group"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;cr1-group;l2-group;l2"
        => [ "LINE", "L2" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;cr1-group;rex-babar"
        => [ "LINE", "Rex-Babar" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;cr1-group;cr1"
        => [ "LINE", "CR1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;cr1-group;cr1;zenon"
        => [ "LINE", "CR1-Zenon" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;r1-group"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;r1-group;r1-subgroup"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;r1-group;r1-subgroup;r1"
        => [ "LINE", "R1" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;r1-group;r1-subgroup;loa"
        => [ "LINE", "R1-LOA" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;r1-group;i-group"
        => [],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;r1-group;i-group;i"
        => [ "LINE", "I" ],
"root;interspersed_repeat;transposable_element;retrotransposed_element;line;group-ii;group-2;r1-like;r1-group;i-group;jockey"
        => [ "LINE", "I-Jockey" ],
    "root;interspersed_repeat;unknown"             => [ "Unknown" ],
    "root;interspersed_repeat;unknown;centromeric" =>
        [ "Unknown", "centromeric" ],
    "root;interspersed_repeat;pseudogene"           => [],
    "root;interspersed_repeat;pseudogene;rna"       => [],
    "root;interspersed_repeat;pseudogene;rna;scrna" => [ "scRNA" ],
    "root;interspersed_repeat;pseudogene;rna;rrna"  => [ "rRNA" ],
    "root;interspersed_repeat;pseudogene;rna;trna"  => [ "tRNA" ],
    "root;interspersed_repeat;pseudogene;rna;snrna" => [ "snRNA" ],
    "root;artifact"                                 => [ "ARTEFACT" ],
    "root;segmental_duplication"                    => [ "Segmental" ],
    "root;tandem_repeat"                            => [],
    "root;tandem_repeat;simple"                     => [ "Simple_repeat" ],
    "root;tandem_repeat;satellite"                  => [ "Satellite" ],
    "root;tandem_repeat;satellite;centromeric"      =>
        [ "Satellite", "centromeric" ],
    "root;tandem_repeat;satellite;macro"         => [ "Satellite", "macro" ],
    "root;tandem_repeat;satellite;y-chromosomal" =>
        [ "Satellite", "Y-chromosome" ],
    "root;tandem_repeat;satellite;acromeric" => [ "Satellite", "acromeric" ],
    "root;tandem_repeat;satellite;w-chromosomal" =>
        [ "Satellite", "W-chromosome" ],
    "root;tandem_repeat;satellite;subtelomeric" =>
        [ "Satellite", "subtelomeric" ],
    "root;low_complexity" => [ "Low_complexity" ]
  };
  my $type = "Unk";
  $type = $cache->{ lc( $class ) }->[ 0 ]
      if ( exists $cache->{ lc( $class ) } );
  my $subtype = "Unk";
  $subtype = $cache->{ lc( $class ) }->[ 1 ]
      if ( exists $cache->{ lc( $class ) } );
  return ( $type, $subtype );
}

1;
