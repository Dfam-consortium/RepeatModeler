#!/usr/local/bin/perl 
##---------------------------------------------------------------------------##
##  File:
##      @(#) rmblast.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A simpler interface to rmblastn
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

=head1 NAME

 rmblast.pl - A simple interface to rmblastn

=head1 SYNOPSIS

  rmblast.pl [-version] [-a(lignments)|-o(outformat)] [-m(atrix=<file>]
             [-gap_init|-gi=<value>] [-gap_ext|-ge=<value>]
             [-minmatch|-mm=<value>] [-minscore|-ms=<value>]
             [-masklevel|-ml=<value>] 
             [-bandwidth|-bw=<value>] [-threads=<value>]
             [-c(leanup] 
             <query_file> <subject_file>

=head1 DESCRIPTION

A simple wrapper for rmblastn which provides sensible defaults for TE
sequence analysis, cross_match-like option names/output format, and finally 
automatic handling of rmblastn database generation.

The options are:

=over 4

=item -version

Displays the version of the program

=item -alignments

Produce alignments for all matches.

=item -outformat

Use the RepeatMasker *.out format.  This only differs from normal cross_match output
in that the orientation field is either "+"/"C" instead of the usual ""/"C".

=item -matrix

Specify the matrix file to use.  Default: ctools.matrix

=item -gap_init

Specify the gap open parameter.  Default: -25

=item -gap_ext

Specify the gap extension parameter.  Default: -5

=item -minmatch

Specify the minimum match parameter.

=item -minscore

Specify the minimum score parameter.  Default: 200

=item -masklevel

Specify the masklevel parameter.  Default: 80

=item -bandwidth

Speicify the bandwidth parameter.

=item -threads

Specify the number of threads to use. Default: 6

=item -cleanup

Remove the database files generated for the subject file.  Otherwise
they remain and are used by subsequent invocations of the script.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2019-2021 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use Getopt::Long;
use POSIX qw(:sys_wait_h ceil floor);
use File::Copy;
use File::Spec;
use File::Path;
use File::Basename;
use Cwd qw(abs_path getcwd cwd);
use Data::Dumper;
use Time::HiRes qw( gettimeofday tv_interval);

# RepeatMasker Libraries
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use SearchResult;
use SearchResultCollection;
use NCBIBlastSearchEngine;

my $Version = $RepModelConfig::VERSION;
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my $version;
my $alignments;
my $matrix = "$FindBin::RealBin/../Matrices/ncbi/nt/ctools.matrix";
my $gapInit = -25;
my $gapExt;
my $insGapExt = -5;
my $delGapExt = -5,
my $minmatch = 7;
my $minscore = 200;
my $masklevel = 80;
my $bandwidth;
my $cleanup;
my $threads = 6;
my $outformat;

my %getopt_args = (
    'version' => \$version, # print out the version and exit
    'alignments|a' => \$alignments,
    'matrix|m=s' => \$matrix,
    'gap_init|gi=i' => \$gapInit,
    'gap_ext|ge=i' => \$gapExt,
    'ins_gap_ext|ige=i' => \$insGapExt,
    'del_gap_ext|dge=i' => \$delGapExt,
    'minmatch|mm=i' => \$minmatch,
    'minscore|ms=i' => \$minscore,
    'masklevel|ml=i' => \$masklevel,
    'bandwidth|bw=i' => \$bandwidth,
    'threads=i' => \$threads,
    'outformat|o' => \$outformat,
    'cleanup|c' => \$cleanup,
);

Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions( %getopt_args)) {
    usage();
}

if ( $outformat && $alignments ) {
  print "\n\nOnly one of '-alignments', '-outformat' may be specified\n\n";
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $version ) {
  print "$Version\n";
  exit;
}

#
# ARGV Processing
#
if ( !@ARGV  ) {
  usage();
}

my $config           = $RepModelConfig::configuration;
my $searchEngineN = NCBIBlastSearchEngine->new(
                       pathToEngine => $config->{'RMBLAST_DIR'}->{'value'} . "/rmblastn" );

$searchEngineN->setGenerateAlignments( 1 ) if ( $alignments );
$searchEngineN->setMatrix( $matrix ) if ( $matrix );
$searchEngineN->setBandwidth( -(abs($bandwidth)) ) if ( $bandwidth );
$searchEngineN->setMaskLevel( $masklevel ) if ( $masklevel );
$searchEngineN->setMinScore( $minscore ) if ( $minscore );
$searchEngineN->setGapInit( $gapInit ) if ( $gapInit );
if ( $gapExt )
{
  $searchEngineN->setInsGapExt( $gapExt );
  $searchEngineN->setDelGapExt( $gapExt );
}
$searchEngineN->setInsGapExt( $insGapExt ) if ( $insGapExt );
$searchEngineN->setDelGapExt( $delGapExt ) if ( $delGapExt );
$searchEngineN->setMinMatch( $minmatch ) if ( $minmatch );
$searchEngineN->setScoreMode( SearchEngineI::complexityAdjustedScoreMode );
$searchEngineN->setCores( $threads ); 

if ( ! -s $ARGV[0] )
{
  die "$ARGV[0] missing or empty!\n";
}

my $query = $ARGV[0];
my $subj = $ARGV[1];

if ( ! -s "$subj.nhr" ) {
  system(   $config->{'RMBLAST_DIR'}->{'value'} . "/makeblastdb -out $subj "
          . "-parse_seqids -dbtype nucl -in $subj > "
          . "makeblastdb.log 2>&1" );
  system("cat makeblastdb.log");
  unlink("makeblastdb.log");
}

$searchEngineN->setQuery( $query );
$searchEngineN->setSubject( $subj );
print STDERR "Search Parameters: " . $searchEngineN->getParameters() . "\n";
my ( $status, $resultCollection ) = $searchEngineN->search();

if ( $status )
{
  print STDERR "\nERROR from search engine (", $? >> 8, ") \n";
} else
{
  for ( my $k = 0 ; $k < $resultCollection->size() ; $k++ ) {
    my $resultRef = $resultCollection->get( $k );
    if ( $alignments ) {
      print "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
    }elsif ( $outformat ) {
      print "" . $resultRef->toStringFormatted( SearchResult::OutFileFormat );
    }else {
      print "" . $resultRef->toStringFormatted();
    }
  }
}

if ( $cleanup ) {
  foreach my $suffix ( 'ndb', 'nos', 'ntf', 'nto', 'not', 'nhr', 'nin', 'nog', 'nsq' ) {
    my $file = $subj . "." . $suffix;
    unlink ($file) if ( -e $file );
  }
}

1;
