#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) fasta-trf-filter.pl
##  Author:
##      Jeb Rosen          jrosen@systemsbiology.org
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Based on Jeb's fasta-trf-filter.py python script.
##

=head1 NAME

fasta-trf-filter.pl - 
    Analyzes the sequences in FASTA file, and prints to standard output only
    the ones that pass the filter.

=head1 SYNOPSIS

  fasta-trf-filter.pl [-min_unmasked_perc # -min_unmasked_contig_bp #]
                      [-max_period #]
                      <fasta_file>

=head1 DESCRIPTION

  Filter sequences from the input file which do not pass the filter after
  being masked by TRF.  TRF is run with the default parameters:

     Match = 2
     Mismatch = 7
     Delta = 7
     PM = 80
     PI = 10
     Minscore = 50
     MaxPeriod = 500
     
  The sequences are then evaluated to determine how much of the input sequence
  remains unmasked and of that how long is the longest unmatched region of the
  sequence.

The options are:

=over 4

=item -min_unmasked_fraction #

The minimum fraction of sequence that must remain unmasked.  Default: 0.2

=item -min_unmasked_contig_bp #

The minimum contiguous length of unmasked sequence. Default: 100

=item -max_period #

The maximum period TRF will look for in the sequence. Default: 500

=back

=head1 SEE ALSO

ReapeatModeler

=head1 AUTHOR

Jeb Rosen <jrosen@systemsbiology.org>

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
my $Version = $RepModelConfig::VERSION;

#
# Paths
#
my $config           = $RepModelConfig::configuration;
my $TRF_PRGM        = $config->{'TRF_PRGM'}->{'value'};

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-version',
                    '-min_unmasked_fraction=s',
                    '-min_unmasked_contig_bp=s',
                    '-max_period=s',
);
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


my $minUnmaskedFractionThreshold = 0.2;
$minUnmaskedFractionThreshold = $options{'min_unmasked_fraction'} if ( $options{'min_unmasked_fraction'} );

my $minUnmaskedContigBPThreshold = 100;
$minUnmaskedContigBPThreshold =  $options{'min_unmasked_contig_bp'} if ( $options{'min_unmasked_contig_bp'} );

my $maxPeriod = 500;
$maxPeriod = $options{'max_period'} if ( $options{'max_period'} );


if ( ! @ARGV || ! -s $ARGV[0] ) {
  usage();
}

my $inputFile = $ARGV[0];

# Output directory
my $outdir = cwd;


my $outFile = "$inputFile.2.7.7.80.10.50.$maxPeriod";

# -ngs : produce compact *.dat output on multisequence files
# -h   : don't make a .html report
my $trfCmd = "$TRF_PRGM $inputFile 2 7 7 80 10 50 $maxPeriod -h -m > $outFile >/dev/null 2>&1";
system($trfCmd);

my %filteredSeqs = ();
if ( -s $outFile ) {
  open IN,"<$outFile.mask" or die "Could not open $outFile.mask for reading!\n";
  my $id;
  my $seq;
  while ( <IN> ) {
    if ( /^>(\S+)/ ) {
      my $tmpID = $1;
      if ( $seq ) {
        my ( $seqLen, $nCount, $maxUnmaskedSubseqLen ) = getStats($seq);          
        my $unmaskedRatio = ($seqLen - $nCount) / $seqLen;
        # If both the sequence is mostly masked and there are no long regions left, filter it!
        if ( $unmaskedRatio < $minUnmaskedFractionThreshold && 
             $maxUnmaskedSubseqLen < $minUnmaskedContigBPThreshold ) {
          $filteredSeqs{$id} = 1;
        }
      }
      $id = $tmpID;
      $seq = "";
      next;
    }
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  if ( $seq ) {
    my ( $seqLen, $nCount, $maxUnmaskedSubseqLen ) = getStats($seq);          
    my $unmaskedRatio = ($seqLen - $nCount) / $seqLen;
    # If both the sequence is mostly masked and there are no long regions left, filter it!
    if ( $unmaskedRatio < $minUnmaskedFractionThreshold && 
         $maxUnmaskedSubseqLen < $minUnmaskedContigBPThreshold ) {
      $filteredSeqs{$id} = 1;
    }
  }
  close IN;
}
unlink("$outFile") if ( -e "$outFile" );
unlink("$outFile.mask") if ( -e "$outFile.mask" );
unlink("$outFile.dat") if ( -e "$outFile.dat" );
unlink("$outFile.summary.html") if ( -e "$outFile.summary.html" );

open IN,"<$inputFile" or die "Could not open $inputFile for reading!\n";
my $id = "";
my $seq = "";
while ( <IN> ) {
  if ( /^>(\S+)/ ) {
    my $tmpID = $1;
    if ( $seq ) {
      unless ( $filteredSeqs{$id} ) {
        $seq =~ s/(\S{50})/$1\n/g;
        $seq .= "\n"
            unless ( $seq =~ /.*\n+$/s );
        print ">$id\n$seq";
      }
    }
    $id = $tmpID;
    $seq = "";
    next;
  }
  s/[\n\r\s]+//g;
  $seq .= $_;
}
if ( $seq ) {
  unless ( $filteredSeqs{$id} ) {
    $seq =~ s/(\S{50})/$1\n/g;
    $seq .= "\n"
       unless ( $seq =~ /.*\n+$/s );
    print ">$id\n$seq";
  }
}
close IN;
 

sub getStats {
  my $seq = shift;

  my $nCount = ($seq =~ tr/Nn/Nn/);
  
  my @wellDefinedBlocks = sort { length($b) <=> length($a) } split(/[Nn]+/,$seq);

  my $maxUnmaskedSubseqLen = 0;
  if ( @wellDefinedBlocks ) {
    $maxUnmaskedSubseqLen = length($wellDefinedBlocks[0]);
  }
  
  return ( length($seq), $nCount, $maxUnmaskedSubseqLen );
}
  
1;
