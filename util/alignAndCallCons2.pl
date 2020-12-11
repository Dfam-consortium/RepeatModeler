#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) alignAndCallCons2.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Arian Smit         asmit@systemsbiology.org
##  Description:
##      A script to generate a MSA from a set of TE copies and call
##      a new consensus.
##

=head1 NAME

alignAndCallCons2.pl - Align TE instances and call consensus

=head1 SYNOPSIS

  alignAndCallCons2.pl -e(lements) <*.fa> | -de(faults) 
                      [-help]

  Example:

    ./alignAndCallCons2.pl -de

        Assumes that a file named "repseq" exists and contains
        all the sequences to be aligned.  This is a deprecated option
        will eventually be removed in favor of the -e parameter.

  or

    ./alignAndCallCons2.pl -e myinstances.fa

=head1 DESCRIPTION

The options are:

=over 4

=item -e(lements) <*.fa>

A multi-sequence FASTA file containing sequences believed to be elements of 
a single TE family.  

=back

=head1 SEE ALSO

ReapeatModeler, dothemsimple.pl, extendcons.pl, createxmoutstk.pl etc..

=head1 AUTHOR

Arian Smit <asmit@systemsbiology.org>

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use Data::Dumper;
use Getopt::Long;
use Cwd;
use File::Spec;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use NCBIBlastSearchEngine;
use CrossmatchSearchEngine;
use SearchResult;
use SearchResultCollection;

# Program version
my $Version = 0.1;

#
# Paths
#
my $phrapDir = "/usr/local/phrap";
my $matrixDir = "$FindBin::RealBin/../Matrices";
my $RMBLAST_DIR = $RepModelConfig::configuration->{'RMBLAST_DIR'}->{'value'};
my $mafftDir = "/usr/local/mafft/bin";
my $defaultEngine = "crossmatch";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-elements|e=s',
                    '-defaults|de',
                    '-help' );
my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0 | less";
  exit( 1 );
}

usage() unless ( ! $options{'help'} && (
                 $options{'defaults'} || ( $options{'consensus'} && $options{'elements'} )));


# Output directory
my $outdir = cwd;

# Override all other methods at picking an outdir if this option is provided.
$outdir = $options{'outdir'} if ( $options{'outdir'} );
# Assumption below is that the directory doesn't already include trailing "/"
$outdir =~ s/\/$//g;
if ( ! -d $outdir ) {
  die "\n\nOutput directory $outdir doesn't exist!\n\n";
}

my $elesFile = "repseq";
if ( exists $options{'elements'} ) {
  if ( -s $options{'elements'} ) {
    $elesFile = $options{'elements'}; 
  }else {
    print "\n\nElements file \"$options{'elements'}\" is missing or empty!\n\n";
    exit(1);
  }
}elsif ( ! -s $elesFile )
{
  print "\n\nDefault elements file \"repseq\" is missing or empty! See -e option for details.\n\n";
  exit(1);
}


my $cmd = "$mafftDir/mafft --localpair --maxiterate 1000 $elesFile > $outdir/mafft.out";
system($cmd);

$cmd = "$FindBin::RealBin/Linup -name consensus $outdir/mafft.out > $outdir/ali";
system($cmd);
$cmd = "$FindBin::RealBin/Linup -name consensus -consensus $outdir/mafft.out > $outdir/rep";
system($cmd);

1;
