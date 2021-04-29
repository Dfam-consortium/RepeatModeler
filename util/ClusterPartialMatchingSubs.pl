#!/usr/bin/perl -w
###---------------------------------------------------------------------------##
###  File:
###      @(#) ClusterPartialMatchingSubs.pl
###  Author:
###      Arian Smit asmit@systemsbiology.org
###  Description:
###      Cluster young subfamilies
###
##******************************************************************************
##*  This software is provided ``AS IS'' and any express or implied            *
##*  warranties, including, but not limited to, the implied warranties of      *
##*  merchantability and fitness for a particular purpose, are disclaimed.     *
##*  In no event shall the authors or the Institute for Systems Biology        *
##*  liable for any direct, indirect, incidental, special, exemplary, or       *
##*  consequential damages (including, but not limited to, procurement of      *
##*  substitute goods or services; loss of use, data, or profits; or           *
##*  business interruption) however caused and on any theory of liability,     *
##*  whether in contract, strict liability, or tort (including negligence      *
##*  or otherwise) arising in any way out of the use of this software, even    *
##*  if advised of the possibility of such damage.                             *
##*                                                                            *
##******************************************************************************

=head1 NAME

 ClusterPartialMatchingSubs.pl - Cluster young partially matching subfamilies

=head1 SYNOPSIS

 ClusterPartialMatchingSubs.pl [-version] [-h(elp)] 
                               <fasta file with related TE copies>

=head1 DESCRIPTION

The script runs cd-hit on a fastafile with sequences of related TE copies and creates consensus sequences 
for each cluster with more than a given number of sequences (default 3). It is efficient in separating 
partially matching subfamilies, when copies are relatively young and entries are restricted to the extent
of the TE copies. Following this script, alignAndCallConsensus.pl should be run with the combined consensi,
both to polish the results but also because sequences have not always been assigned to the correct cluster
and other sequences may not have been included at all.

The options are:

=over 4

=item -version

Displays the version of the program

=item -n(umber) #   

Minimum number of sequence in a cluster to build a consensus for (3 or higher only makes sense)

=item -lmax #       

Maximum length of a query sequence to include (e.g. given a mix of solo LTRs and full ERVs)

=item -lmin #       

Minimum length of a query sequence to include

=item -de(bug)      

Files will be saved that are by default removed. Among these are stderr text and intermediate 
files by all routines and alignment files created by cd-hit.

=back

cd-hit options: (see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide for details)

=item -mi(nid) #    

Lowest sequence identity (in %) of copies with a representative to cluster with it.

=item -ma(xid) #    

The script performs hierarchical clustering, starting with "maxid" (default 90%) 

=item -s(tep) #     

As lowest identity, lowering a "step" (def. 10%) each round until "minid" (default 80%) 
Expect ave 10% subst TE copies to be 10-25% similar to a subfamily representative.
Too high a minid will give many unlinked seqs, too low will give fail to separate subfams.
Setting -maxid and -minid to the same level performes just a single clustering.

=item -b(andwidth) # 

Default 20. Higher bandwidths could possibly lump together more subs with large indel diffs.

=item -f(ast)       

(-g 0 in cd-hit) Faster cd-hit algorithm. Perhaps useful when optimizing parameters on large sets.

=item -t(hreads) #

Default 4; with 0 all CPUs will be used.

=back

consensus-building options: 

A matrix optimal for 14% subst level from original in a 41% GC background is always used.
We think cd-cluster may perform poorly finding subs for copies that are very much older.
Otherwise default parameters for alignAndCallConsensus.pl

=item -a(uto_off)   

Skips the one AutoRunBlocker run for each cluster (the slow step, especially for long sequences)

=item -di(vergence) % 

By default a substitution matrix and gap penalties are used appropriate for the alignment 
of neutrally evolved copies to the ancestral state with a substitution level of 14%. Higher 
substitution level parameters can be chosen.

=item -g(clevel) %  

By default a substitution matrix optimal for a human-average 41% GC-background is used.

=item -p(rune) #    

By default the MSA will be pruned at both ends until the higher of seqs_in_cluster/3 or 3 
sequences are aligned at a position. -p 4 will add 4 to this number, -p -1 will subtract 1 (never <2)

=item -w(indow) #   

AutoRunBlocker.pl window size. Default 10 (smaller numbers faster and better on older elements) 

=item -z (b#/l#/f#) 

*** the z option NOT yet implemented ***
After the first round of consensus calling and pruning a string of # (def. 7) Zs/Hs is added
to the ends to allow extension in the next rounds. By default no pads are added. 
-z b adds 7 to both sides, -b l10 10 to the left only, -b r12 12 to the right only.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2019-2021 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Arian Smit <asmit@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Storable qw(nstore retrieve);
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
#
use RepModelConfig;
use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
use NCBIBlastSearchEngine;
use SearchResult;
use SearchResultCollection;

# Program version
my $Version = $RepModelConfig::VERSION;

#
# Paths
#
my $config           = $RepModelConfig::configuration;
my $CDHIT_PRGM       = $config->{'CDHIT_DIR'}->{'value'} . "/cd-hit";
my $CDHIT_CLSTR_PRGM = $config->{'CDHIT_DIR'}->{'value'} . "/clstr_rev.pl";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-version',
                     '-auto_off|a',
                     '-bandwidth|b=i',
                     '-debug|de',
                     '-divergence|di=i',
                     '-fast|f',
                     '-gclevel|g=i',
                     '-lmax=i',
                     '-help|h',
                     '-lmin=i',
                     '-maxid|ma=i',
                     '-minid|mi=i',
                     '-number|n=i',
                     '-prune|p=i',
                     '-step=i',
                     '-threads|t=i',
                     '-window|w=i', 
                     '-z=s',
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

if ( exists $options{'help'} ) {
  usage();
}

if ( ! @ARGV || ! -s $ARGV[0] ) {
  die "\n\nFasta file missing!\nUse '$0 -h' to view the help.\n\n";
}

my $fafile = $ARGV[0];

my $number = 3; # default minimum number of sequences in a cluster to create a consensus for.
if ($options{'number'}) {
  if ($options{'number'} < 3) {
    print STDERR "Minimum number to make a consensus for rest to 3. For smaller clusters see the .clstr output file.\n";
    $number = 3;
  } else {
    $number = $options{'number'};
  }
}

my $maxid = 0.90;
$maxid = $options{'maxid'}/100 if $options{'maxid'};
my $minid = 0.80;
$minid = $options{'minid'}/100 if $options{'minid'};
my $step = 0.1;
$step = $options{'step'}/100 if $options{'step'};

my $band = ""; # default of cd-hit is 20
$band = $options{'bandwidth'} if $options{'bandwidth'};

my $accurate = "-g 1";
$accurate = "" if $options{'fast'};

# currently fixed in command line:
# -d 0 : no limit to length of sequence name; default is 20

# -M 1600 : memory limit to 1.6 GB (default is 800 MB) Probably more than enough, 
# but could be set to unlimited using -M 0

my $threads = 4;
$threads = $options{'Threads'} if $options{'Threads'};

# I think -l "length of throw_away_sequences" (default 10) refers to sequence entries that are ignored 
# we could set that considerably higher, as anything < 30 bp or even < 50 bp does not contribute here

# -t "tolerance for redundance" (default 2) must refer to identical entries. We elminate these with 
# DeleteAlmostDupsfromSeqfile.pl beforehand anyway, so probably does not need to be set

# When using -sc 1 the output clusters by decreasing size instad of by
# length of representative (?) and the results were poorer for my
# first test, so I've left that default

# -n wordlength: interesting, but cd-hit does not allow going above word_length 5!
# I was going to lower the wordlength by divergence from 8 down, but 5 is small enough for anything

# alignment parameters, i.e. picking a matrix (alignAndCallConsensus.pl adjusts the gap penalties)

my $matrix = "14";
if ($options{'divergence'} && $options{'divergence'} > 15) {
  $matrix = 18;
  # gap penalty adjustements are done in alignAndCallConsensus.pl
  if ($options{'divergence'} > 18) {
    $matrix = 20;
    if ($options{'divergence'} > 22) {
      $matrix = 25;
    }
  }
}
if ($options{'gclevel'}) {
  my $tempmatrix = $matrix;
  if ($options{'gclevel'} < 40) {
    $tempmatrix = "$matrix"."p39";
    if ($options{'gclevel'} < 38) {
      $tempmatrix= "$matrix"."p37";
      if ($options{'gclevel'} < 36) {
	$tempmatrix= "$matrix"."p35";
      } 
    }
  } elsif ($options{'gclevel'} > 42) {
    $tempmatrix= "$matrix"."p43";
    if ($options{'gclevel'} >44) {
      $tempmatrix= "$matrix"."p45";
      if ($options{'gclevel'} > 46) {
	$tempmatrix= "$matrix"."p47";
	if ($options{'gclevel'} >48) {
	  $tempmatrix= "$matrix"."p49";
	  if ($options{'gclevel'} > 50) {
	    $tempmatrix= "$matrix"."p51";
	    if ($options{'gclevel'} >= 53) {
	      $tempmatrix= "$matrix"."p53";
	    }
	  }
	}
      }
    }
  }
  $matrix = $tempmatrix;
}

my $epoc = time();

my %seq;
my $sqname;
my @seqnames;
open (IN, "$fafile");
while (<IN>) {
  if (/^\s*>(\S+)/) {
    $sqname = $1;
    $seq{$sqname} = $_;
    push @seqnames, $sqname;
  } elsif ($sqname) {
    $seq{$sqname} .= $_ if /\S/;
    # following not necessary if cd-hit has a careful sequence checker
    die "This line in the file $fafile contains (a) non-DNA character(s):\n$_" if /[^acgtnrykmswbdhvzACGTNRYKMSWBDHVZ\s]/;
  } else {
    die "Something is awry before or with this line in the file $fafile\n$_\n";
  }
}   
if ($options{'lmax'} || $options{'lmin'}) {
  my $faselection = "$fafile.selection$epoc";
  my $survivorcount =0;
  open (OUTF, ">$faselection") || die "Cannot write to $faselection";
  foreach my $seqname (@seqnames) {
    my $string = $seq{$seqname};
    $string =~ s/^>.+\n//;
    $string =~ s/\s//g;
    my $len = length $string;
    if ($options{'lmax'} && (length $string) > $options{'lmax'}) {
      print STDERR "$seqname ($len bp) excluded from cluster analysis as it is longer than $options{'lmax'} bp.\n";
    } elsif ($options{'lmin'} && (length $string) < $options{'lmin'}) {
      print STDERR "$seqname ($len bp) excluded from cluster analysis as it is shorter than $options{'lmin'} bp.\n";
    } else {
      print OUTF "$seq{$seqname}";
      ++$survivorcount;
    }
  }
  die "No sequences survived the size selection indicated by -lmin and/or lmax\n" unless $survivorcount;
  $fafile = $faselection if $#seqnames >= $survivorcount;
}


if ($options{'debug'}) {
  open (OUTCL, ">commandlinesdone");
}
my $id = $maxid;
my $infile = $fafile;
my $combinedcluster;
my $round = 0;
while ($id >= $minid) {
  my $tag = 100*$id;
  my $outfile = "$fafile\_$tag";
  print OUTCL "$CDHIT_PRGM -i $infile -o $outfile -c $id $band -M 1600 -T $threads -d 0 $accurate\n" if $options{'debug'};
  system "$CDHIT_PRGM -i $infile -o $outfile -c $id $band -M 1600 -T $threads -d 0 $accurate >> cdhitnoises$epoc";
  if ($round) {
    my $firstclstr = "$infile.clstr";
    if ($combinedcluster) {
      $firstclstr = $combinedcluster;
      $combinedcluster =~ s/\.clstr$/-$tag\.clstr/;
    } else {
      $combinedcluster = "$infile"."-$tag.clstr";
    }
    die "firstclstr ( $firstclstr ) does not exist\n" unless -s $firstclstr;
    die "outfile.clstr ( $outfile.clstr ) does not exist\n" unless -s "$outfile.clstr";
    print OUTCL "$CDHIT_CLSTR_PRGM $firstclstr $outfile.clstr > $combinedcluster\n" if $options{'debug'};
    system "$CDHIT_CLSTR_PRGM $firstclstr $outfile.clstr > $combinedcluster";
    unlink "$firstclstr","$outfile.clstr" unless $options{'debug'};
  }
  my $diff = $id - $minid;
  if ($diff) {
    if ($diff > $step) {
      $id -= $step;
    } else {
      $id = $minid;
    }
  } else {
    last;
  }
  $infile = $outfile;
  ++$round;
}
unlink "cdhitnoises$epoc" unless $options{'debug'};

my ($cluster,$count,$representative) = ();
my $curdir = `pwd`;
chomp $curdir;
if ($maxid > $minid) {
  open (IN, $combinedcluster);
} else {
  my $tag = 100*$maxid;
  my $clusterfile = "$fafile\_$tag.clstr";
  open (IN, $clusterfile);
}

if (-s "clusterconsensi") {
  rename "clusterconsensi","consensi_before$epoc";
}
open (OUTCONS, ">clusterconsensi");
while (<IN>) {
  if (/^>Cluster\s+(\d+)/) {
    my $nr = $1;
    if ($count && $count >= ($number-1)) { #count is zero-based, $number human-based
      &makeconsensus;
    }
    @seqnames = ();
    $cluster = "Cluster$nr"; 
  } elsif (/^(\d+)\s.+>(\S+)\.\.\.\s+(\S+)/) {
    # I think those three dots always finish the sequence name in cd-hit .clstr files
    $count = $1;
    push @seqnames, $2;
    $representative = $2 if $3 eq '*';
  }
}
if ($count && $count >= ($number-1)) {
  &makeconsensus;
}
close OUTCONS;

sub makeconsensus {
  print STDERR "Creating consensus for $cluster\n";
  system "mkdir $cluster";
  chdir $cluster;
  open (OUT, ">repseq");
  foreach my $seqname (@seqnames) {
    print OUT "$seq{$seqname}";
  }
  close OUT;
  my $seqnr = $#seqnames + 1;
  open (OUT, ">rep");
  $seq{$representative} =~ s/^>/>$cluster ($seqnr sequences) /;
  print OUT "$seq{$representative}";
  close OUT;
  my $minaligned = int($seqnr/3);
  $minaligned = int($seqnr/2.5) if $minaligned < 3;
  $minaligned = 2 if $minaligned < 2;
  if ($options{'prune'}) {
    $minaligned += $options{'prune'};
  }
  $minaligned = 2 if $minaligned < 2;
  print OUTCL "$cluster\]\$ $FindBin::RealBin/alignAndCallConsensus.pl -ma $matrix -p $minaligned -re 3 -rm\n" if $options{'debug'};
  system "$FindBin::RealBin/alignAndCallConsensus.pl -de -ma $matrix -p $minaligned -re 3 -rm >consensusbuildingnoises1 2>&1";
  if ($options{'auto_off'}) {
    system "$FindBin::RealBin/alignAndCallConsensus.pl -de -ma $matrix -p $minaligned -rm -ht >consensusbuildingnoises 2>&1";
  } else {
    my $ws = 10;
    $ws = $options{'window'} if $options{'window'};
    system "$FindBin::RealBin/AutoRunBlocker.pl -linup ali -windowSize $ws -minRatioAgreement 1.5 -minCopyAgreement $minaligned >tempout 2>&1";
    my $tempnew = "temp$epoc";
    open (INTEMP, "tempout");
    my $print = "";
    while(<INTEMP>) {
      if (/There are no changes/) {
	last;
      } elsif (/^>/) {
	open (OUT, ">$tempnew");
	print OUT ">$cluster $representative\n";
	$print = 1;
      } elsif ($print && /\w/) {
	print OUT;
      }
    }
    close OUT;
    close INTEMP;
    if (-s $tempnew) {
      system "mv $tempnew rep";
      system "$FindBin::RealBin/alignAndCallConsensus.pl -de -ma $matrix -re 2 -rm -ht";
    } else {
      system "$FindBin::RealBin/alignAndCallConsensus.pl -de -ma $matrix -rm -ht";
    }
  }
  unlink "repseq.log","tempxmatch.stderr", "tempout", "repbefore1" unless $options{'debug'};
  open (INREP, "rep");  
  while (<INREP>) {
    print OUTCONS if /\S/; 
    # somehow alignAndCallConsensus.pl prints out empty lines before and
    # after the sequence; part of a bug that Robert is fixing anyhow
  }
  close INREP;
  chdir $curdir;
}


