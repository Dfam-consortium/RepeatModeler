#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) TSD.pl
##  Author:
##      Arian Smit         asmit@systemsbiology.org
##  Description:
##      Identify putatative TSD lengths in a TE family
##

=head1 NAME

TSD.pl - Identify putative TSD lengths in a TE family

=head1 SYNOPSIS

  TSD.pl <linup *.ali file>

=head1 DESCRIPTION

 input file is a list of "name position "flanking-bases" derived from an extracttomultaln or Linup .ali output file, like:
 NOTE: At some point the ali files were reworked to have only 1 space after the index rather than 2.  

    1_152_107531_108     187 GATGGGCC
    1_22_305123_3072     174 AGCCCAGA
    11_44_6264_8453       96 TAGGTCTC
    12_10_964177_965     186 TGCTGGAT
    ...

    1_152_107531_108     913 GGCCCTAA
    11_44_6264_8453     1964 TCCAAGCA
    12_10_964177_965     840 TGGGCATG
    12_23_157980_159     797 CAAACACA

 copies don't have to have a match extending to both beginning and end
 sorting by name is not necessary (though can't hurt unless somehow the end ends up prceding the beginning of a copy)
 length of flanking bases on either side should be the same

The options are:

=over 4

=item <linup *.ali file>

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
my $Version = $RepModelConfig::VERSION;

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0 | less";
  exit( 1 );
}

usage() unless ( @ARGV );

my ($nextisleft, $nextisright, $extend, $skip) = ();
my $lengthPad = 0;
my $lengthPadend = 0;
my %beg = ();
my %end = ();
while (<>) {
  my @bit = split;
  if ($bit[2]) {
    if ($bit[2] =~ /^([ZH]+)[ACGTNX]/ && !$extend) {
      $lengthPad = length $1;
      $nextisleft = 1;
    } elsif ($nextisleft) {
      if (/\d\s([ACGTNX]{$lengthPad})/) {
	$beg{$bit[0]} = $1;
      } else {
	# too many spaces after digit means that copy did not extend to the beginning
	$nextisleft = 0;
      }
      #other end
    } elsif ($bit[2] =~ /(\S*?)([ZH]+)$/ && !$extend) {
      $skip = length $1;
      $lengthPadend = length $2;
      $extend = 1 if $lengthPadend < $lengthPad;
      $nextisright = 1;
    } elsif ($bit[2] eq /([ZH]+){$lengthPad}/ && !$extend) {
      $skip = 0;
      $lengthPadend = length $2;
      $nextisright = 1;
    } elsif ($nextisright) {
      if (/\d \S{$skip}([ACGT]{$lengthPadend})/) {
	$end{$bit[0]} = $1;
      }
    } elsif ($extend) {
      $end{$bit[0]} .= $bit[2] if $end{$bit[0]};
    }
  } else {
    $nextisleft = $nextisright = $skip = 0;
  }
}

my %score = ();
foreach my $name (sort keys %beg) {
  next unless $end{$name};
  my @left = split "", $beg{$name};
  my @right = split "", $end{$name};
      print "$beg{$name} <=> $end{$name}\n";
  for (my $i = 0; $i <= $#left; ++$i) {
    my $length = $i+1;
    for (my $j = 0; $j <= $i; ++$j) {
      my $k = $#left-$i+$j;
      $score{$length} += &calc($left[$k], $right[$j]) if $left[$k] && $right[$j];
    }
    # given a 42% GC level (though it's unclear if target sites for DNA
    # transposons and LTR elements are more or less AT-rich than the
    # rest of the genome), the average score expected per base is
    # 9*.29-6*.21-15*.21-17*.29 = -6.73 for [AT] and -7.14 for
    # [GC]. If background is 45% GC, the scores are both -6.9.
    $score{$length} += $length*(6.9);
  }
}

foreach my $len (sort keys %score) {
  my $score = int($score{$len}/$len+0.5);
  print "$len $score\n";
}

sub calc {
  my $l = shift;
  my $r = shift;
  if ($l eq "N" || $r eq "N" ){ 
    return 0;
  } elsif ($l eq $r) {
    if ($l =~ /[AT]/) {
      return 9;
    } else {
      return 10;
    }
  } else {
    if ($l =~ /[AG]/ && $r =~ /[AG]/ ||
	$l =~ /[CT]/&& $r =~ /[CT]/) {
      return "-6";
    } elsif ($l =~ /[AT]/ && $r =~ /[AT]/) {
      return "-17";
    } else {
      return "-15";
    }
  }
}
