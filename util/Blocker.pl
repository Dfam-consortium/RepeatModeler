#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) Blocker.pl
##  Author:
##      Arian Smit         asmit@systemsbiology.org
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Identify gaps (insertions/deletions) in a consensus
##

=head1 NAME

Blocker.pl - Identify insertions/deletions in a Linup region

=head1 SYNOPSIS

  Blocker.pl -linup <ali file>
             -startBlock #
             -startColumn #
             -endBlock #
             -endColumn #
             [-cmdiffs]

  or

  Blocker.pl <linup_ali_file> <startBlock> <startColumn> <endBlock> <endColumn>

  Examples:

     # Analyse the region in the block starting with consensus position 35 and
     # from the 20th column to the end of the block.
     ./Blocker.pl ali 35 20
   or
     ./Blocker.pl ali 35 20 35

     # Analyse the region in the block starting with consensus position 35 and
     # from the 20th column to the 5th column of the block starting with consensus
     # position 59.
     ./Blocker.pl ali 35 25 59 5
     
     # Analyse the region in the block starting with consensus position 35 and
     # from the 20th column to the end of the block starting with consensus
     # position 59.
     ./Blocker.pl ali 35 25 59 0
   or 
     ./Blocker.pl ali 35 25 59

     # Analyse the region in the block starting with consensus position 59 and
     # from the 3rd column to the 7th column in from the end of the block starting
     # with consensus position 113.
     ./Blocker.pl ali 59 3 113 -7

     # Same as above but using named parameters
     ./Blocker -linup ali -startBlock 59 -startColumn 3 -endBlock 113 -endColumn -7 

=head1 DESCRIPTION

This is a specialized tool for considering regions of a Linup output file 
( aka "ali" file ).  A typical Linup file is formatted into screen formatted
portions of a larger alignment.  Each section or "block" looks like:

  consensus         1 TG----G----AC----T----G 7
  ref:seq1          1 TG----G----AC----T----G 7
  seq2             10 TA----G----AC---------G 16
  seq3             12 TG----G----AC----T----G 18
  seq4             15       G----AC----TA---G 20 

A block is terminated with a blank line and the next "block" continues with 
the next section of the alignment. For example:

  consensus         8 --G-G-CGTA-T----GCC---A 18
  ref:seq1          8 A-G-G-C-TA-T----GCC---A 18
  seq2             17 --G-G-GGTA-TACGTGCC---A 31
  seq3             19 --G-G-CGAACT----G-C---A 29
  seq4             21 --A-G-CGTA-T----GCC---A 31

The "startBlock"/"endBlock" parameters refer to the first consensus position
in a Linup block.  To refer to the second block in the example above one 

=head1 DESCRIPTION

This is a specialized tool for considering regions of a Linup output file 
( aka "ali" file ).  A typical Linup file is formatted into screen formatted
portions of a larger alignment.  Each section or "block" looks like:

  consensus         1 TG----G----AC----T----G 7
  ref:seq1          1 TG----G----AC----T----G 7
  seq2             10 TA----G----AC---------G 16
  seq3             12 TG----G----AC----T----G 18
  seq4             15       G----AC----TA---G 20 

A block is terminated with a blank line and the next "block" continues with 
the next section of the alignment. For example:

  consensus         8 --G-G-CGTA-T----GCC---A 18
  ref:seq1          8 A-G-G-C-TA-T----GCC---A 18
  seq2             17 --G-G-GGTA-TACGTGCC---A 31
  seq3             19 --G-G-CGAACT----G-C---A 29
  seq4             21 --A-G-CGTA-T----GCC---A 31

The "startBlock"/"endBlock" parameters refer to the first consensus position
in a Linup block.  To refer to the second block in the example above one 
would use "-startBlock 8".  To refer to both the first and second blocks
one would use "-startBlock 1 -endBlock 8".

The column parameters refer to the column number within a block starting from 1
( consensus or gap ).  For instance the parameters "-startBlock 8 -startColumn 4"
point to the gap column between the two "G" consensus calls in the second
block ( "G-G" ).  In addition, a column index of 0 is mapped to the last column
in a block and negative columns are counted backwards from the end of a block.

note: script goes awry when not all entries have a unique name, including the 
name consensus

The options are:

=over 4

=item -linup <ali file>

blah

=back

=head1 SEE ALSO

ReapeatModeler, extendcons.pl, createxmoutstk.pl, dothemsimple.pl etc.

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
#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};


# Program version
my $Version = $RepModelConfig::VERSION;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-linup=s',
                    '-startBlock=s',
                    '-endBlock=s',
                    '-startColumn=s',
                    '-endColumn=s',
                    '-cmdiffs'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit( 1 );
}

my $alifile;
my $beginblock;
my $begin;
my $endblock;
my $end = 0;
# Two usages....legacy is:
#     ./Blocker <linupFile> <startBlock> <startColumn> [<endBlock>] [endColumn]
# And current is:
#     ./Blocker -linup ali -startBlock 1 -startColumn 10 [-endBlock #] [-endColumn #] 
if ( $options{'linup'} && $options{'startBlock'} && $options{'startColumn'} ) 
{
  $alifile = $options{'linup'}; 
  $beginblock = $options{'startBlock'};
  $begin = $options{'startColumn'};
  $endblock = $beginblock;
  $endblock = $options{'endBlock'} if ( $options{'endBlock'} ); 
  $end = $options{'endColumn'} if ( $options{'endColumn'} );
}elsif ( @ARGV >= 3 ){
  $alifile = shift @ARGV;
  $beginblock = shift @ARGV;
  $begin = shift @ARGV;
  $endblock = $beginblock;
  $endblock = shift @ARGV if ( @ARGV );
  $end = shift @ARGV if ( @ARGV );
}else{
  usage();
}

if ( $begin <= 0 ) {
  print "\n\nstartColumn must be > 0!\n\n";
  usage();
}
if ( $endblock < $beginblock ) {
  print "\n\nendBlock must be >= startBlock!\n\n";
  usage();
}

# For the purposes of labeling substitutions.
my %mutChar = (
                "CT" => 'i',
                "TC" => 'i',
                "AG" => 'i',
                "GA" => 'i',
                "GT" => 'v',
                "TG" => 'v',
                "GC" => 'v',
                "CG" => 'v',
                "CA" => 'v',
                "AC" => 'v',
                "AT" => 'v',
                "TA" => 'v' );


$end = 0 unless $end;
my $blocktoprint = "";
my $span = 0;
my $blockwidth = 0;
my $take = "";
my %string = ();
my $spacelen = 0;
my %done;
my %dontdo;
my $pos;
my $nextline;
my $prevline;
my $prevcons;
open (IN, $alifile) or die;
while (<IN>) {
  last if (/^>/);
  if (/^consensus\s+(\d+)(\s+)(\S+)\s+(\d+)/) {
#    last if /^consensus sequence:/;  Must have been an old format
    my $beg = $1;
    $spacelen = length $2;
    $blockwidth = length $3;
    %done = ();
    next if $beg < $beginblock;
    last if $beg > $endblock;
    $take = 1;
    $pos = 0;
    $span = $blockwidth;
    if ($beg == $beginblock) {
      $pos = $begin - 1;
      $span -= $pos;
    }
    if ($beg == $endblock) {
      $span += $end if $end < 0;
      $span = $end if $end > 0;
      $span = $end - $begin + 1 if ($beginblock == $endblock && $end > 0);
#only fails in pathological case when someone types a negative beginning and a positive end
#they deserve it
    }
  }
  if ( $take && /^(\S+)\s+\d+\s{$spacelen}(.{$blockwidth})/ ) {
    my $name = $1;
    my $string = $2;
    my $cnt = 0;
    while ($done{$name}) {
      ++$cnt;
      $name .= "_$cnt";
    }
    $done{$name} = 1;
    next if $dontdo{$name};
    if ($nextline) {
      $nextline = 0;
      $prevcons = $name;
    }	
    $nextline = 1 if $name eq 'consensus';
    my $sub = substr($string,$pos,$span);
    if ($sub =~ /\s+/) { # alignment starts or ends after block begins or ends
      $dontdo{$name} = 1;
      $string{$name} = ""; # if alignment started at same spot
      next;
    }
    $sub =~ tr/-//d;
    $string{$name} .= $sub;
  }
}

my $conslen = length $string{'consensus'};

my $maxname = 0;
my $maxsub = 0;
my %length;
foreach my $key (keys %string) {
  $length{$key} = (length $string{$key});
  $maxsub = $length{$key} if $length{$key} > $maxsub;
  $maxname = (length $key) if (length $key) > $maxname;
}

my $space = " " x ($maxname - 9) ; # consensus is 9 characters;
$blocktoprint = "consensus $space $string{'consensus'}\n";

$space = " " x ($maxname - (length $prevcons)) ; 
$blocktoprint .= "$prevcons $space $string{$prevcons}\n";
$length{'consensus'} = $length{$prevcons} = 0; 

my $bestlength = 0;
my $maxcnt = 0;
my $conscnt = 0;
my $secondbestcnt = 0;
my $secondbestlength = 0;
my %cnt = ();
my $totalcnt = 0;
for (my $i = 1; $i <= $maxsub; $i++) {
  foreach my $name (keys %length) {
    if ($length{$name} == $i) {
      my $diff = $maxname - (length $name);
      $space = " " x $diff;
      $blocktoprint .= "$name $space $string{$name}\n";
      ++$cnt{$i};
      ++$totalcnt;
    }
  }
  if ($cnt{$i}) {
    $conscnt = $cnt{$i} if $conslen == $i;
    if ($cnt{$i} >= $maxcnt) {
      $secondbestlength = $bestlength;
      $bestlength = $i;
      $secondbestcnt = $maxcnt;
      $maxcnt = $cnt{$i};	
    } elsif ($cnt{$i} > $secondbestcnt) {
      $secondbestcnt = $cnt{$i};
      $secondbestlength = $i;
    }
  }
}

# to report a cluster of longer sequences; these tend to be off the screen
my $outlier = "";
if ($maxsub > $bestlength + 3) {
  my $diffinlen = 0;
  my $diffinlen2 = $secondbestlength - $bestlength;
  for (my $i = $bestlength + 3; $i <= $maxsub; $i++) {
    my $lengthdiff = $i - $bestlength;
    if ($cnt{$i} && $cnt{$i} > 3 && $cnt{$i} > $maxcnt/10 &&
      $i - $secondbestlength > 2 && $lengthdiff > $diffinlen2) {
      $diffinlen = $i - $bestlength;
      $outlier = ", $cnt{$i} copies of $i bp";
    }
  }
}

my %col;
foreach my $name (keys %length) {
  if ($length{$name} == $bestlength) {
    my @string = split "", $string{$name};
    for (my $i = 0; $i < $bestlength; ++$i) {
      $col{$i} .= $string[$i];
    }
  }
}

open (IN, "$FindBin::RealBin/../Matrices/linupmatrix ") or die;
my @cons;
my @l;
my %matrix;
while (<IN>) {
  s/^\s+//;
  if (/^A/) {
    @l = split;
    @cons = @l;
  } elsif (/^[\d-]/) {
    my @s = split;
    my $l = shift @cons;
    for (my $i = 0; $i <= $#s; ++$i) {
      $matrix{$l}{$l[$i]} = $s[$i];
    }
  }
}
close IN;

my $printalert = "";
$printalert = "Length difference. " unless $bestlength == (length $string{'consensus'});
my @oldstring = split "", $string{'consensus'};
my $newstring = "";
my $maxl;
my $spacer = "    ";
for (my $i = 0; $i < $bestlength; ++$i) {
  my @s = split "", $col{$i};
  my $max = -999999999;
  foreach my $let ( 'A','C','G','T','N') {
    my $score = 0;
    for (my $j = 0; $j <= $#s; ++$j) {
      $score += $matrix{$let}{$s[$j]};
    }
    if ($score > $max) {
      $maxl = $let;
      $max = $score;
    }
  }
  $newstring .= "$maxl";
  if ($i <= $#oldstring) {
    if ($maxl ne 'N' && $oldstring[$i] eq 'N' ||
	$maxl =~ /[AG]/ && $oldstring[$i] =~ /[CT]/ ||
	$maxl =~ /[CT]/ && $oldstring[$i] =~ /[AG]/) {
      $printalert = "Sequence difference. " unless $printalert;
      $spacer .= $maxl;
    } else {
      $spacer .= " ";
    }
  }
}

# Calculate a visual diff line between the old/new consensus
my $diff = $spacer;
if ( $options{'cmdiffs'} ) {
  $diff = "    ";
  my $diffLen = length($string{'consensus'});
  $diffLen = length($newstring) if ( length($newstring) > $diffLen );
  for ( my $i = 0; $i < $diffLen; $i++ )
  {
    my $newB;
    $newB = substr( $newstring, $i, 1 ) if ( $i < length($newstring) );
    my $oldB;
    $oldB = substr( $string{'consensus'}, $i, 1 ) if ( $i < length($string{'consensus'}));
    if ( ! defined $newB || ! defined $oldB ) {
      $diff .= "-";
    }elsif ( my $mc = $mutChar{ uc( $newB . $oldB ) } )
    {
      $diff .= $mc;
    }elsif (    ( $newB =~ /[BDHVRYKMSWNX]/i )
             || ( $oldB =~ /[BDHVRYKMSWNX]/i ) )
    {
      $diff .= "?";
    }else
    {
      $diff .= " ";
    }
  }
}

if ($printalert) {
  if ($printalert =~ /^Len/) {
    print "\n$printalert$maxcnt of $totalcnt for $bestlength bp ($secondbestcnt second best for $secondbestlength bp$outlier, $conscnt for original $conslen bp)\nold $string{'consensus'}\n$diff\nnew $newstring\n";
  } else {
    print "\n$printalert$maxcnt of $totalcnt for $bestlength bp ($secondbestcnt second best for $secondbestlength bp$outlier)\nold $string{'consensus'}\n$diff\nnew $newstring\n";
  }
} else {
  print "\n$maxcnt of $totalcnt for $bestlength bp ($secondbestcnt second best for $secondbestlength bp$outlier)\nsame $string{'consensus'}\n";
}
print "\n$blocktoprint\n";
