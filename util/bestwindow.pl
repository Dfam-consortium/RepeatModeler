#!/usr/bin/perl -w

# Finds a long region in a consensus to which many copies align.
# By default it optimizes the product of copy # and window size..
# Change $copiesonly parameter to just get highest copy number (which runs a lot faster)

die "bestwindow.pl <cross_match out file> minimum-windowsize 
Any third argument will make it count copies instead of bases\n" unless $ARGV[1];

my $file = shift;
my $minimumwindow = shift;
my $copiesonly = 0;
$copiesonly = 1 if $ARGV[0];

open (IN, $file);
#cross_match output file batch of queries vs one consensus
my %beg;
my %end;
my $id = 1;
while (<IN>) {
  if (/\(\d+\)\s+\S/) {
    my @bit = split;
    if ( $bit[8] eq 'C' && $bit[10] =~ /^\(/ ) {
      $beg{$id} = $bit[12];
      $end{$id} = $bit[11];
    } elsif ( $bit[8] eq 'C' && $bit[11] =~ /^\(/ ) {
      $beg{$id} = $bit[13];
      $end{$id} = $bit[12];
    } elsif ( $bit[8] eq '+' && $bit[11] =~ /^\d+$/ ) {
      $beg{$id} = $bit[11];
      $end{$id} = $bit[12];
    } elsif ( $bit[9] =~ /^\d+$/ ) {
      $beg{$id} = $bit[9];
      $end{$id} = $bit[10];
    } else {
      die "This is wrong\n$_";
    }
    next if $end{$id} - $beg{$id} <= $minimumwindow;
    # This short match will never be counted (no window is comprised
    # in it) so might as well leave it out. It will be "overwritten" by
    # the next (if that's long enough) as the $id doesn't change now.
    ++$id;
  }
}

my @matches;
my @beg;
my @end;
my %donebeg;
my %doneend;

# Bloodboilingly, the use of &subroutine (&byendbegin) has been
# abandoned or at least is not functioning properly in the following:
foreach (sort byendbegin keys %beg) {
  push @matches, $beg{$_};
  push @matches, $end{$_};
  unless ( $donebeg{$beg{$_}} ) {
    push @unsortedbeg,$beg{$_};
    ++$donebeg{$beg{$_}};
  }
  unless ( $doneend{$end{$_}} ) {
    push @end, $end{$_};
    ++$doneend{$end{$_}};
  }
}

sub byendbegin {
  $end{$a} <=> $end{$b} ||
      $beg{$a} <=> $beg{$b};
}

@beg = sort {$a<=>$b} @unsortedbeg;
# @end is already sorted

my $score = ();
for (my $i = 0; $i <= $#beg; ++$i) {
  print STDERR "$beg[$i] of $beg[$#beg]\n";
  while ($end[0] && $end[0] - $beg[$i] < $minimumwindow - 1) {
    print STDERR "End $end[0] shifted because begin is $beg[$i]\n";
    shift @end;
  }
  last unless $end[0];
  #sorted by end
  while ($matches[1] < $beg[$i] + $minimumwindow - 1) {
    shift @matches;
    shift @matches;
  }
  my @tempmatches = @matches;
  if ($copiesonly) {
    my $window = $beg[$i]."_$end[0]";
    while ($tempmatches[1] < $end[0]) {
      shift @tempmatches;
      shift @tempmatches;
    }
    next if !$tempmatches[0];
    for (my $k = 0; $k <= $#tempmatches; ++$k) {
      ++$score{$window} if $tempmatches[$k] <= $beg[$i];
    }
  } else {
    for ($j = 0; $j <= $#end; ++$j) {
#    print "$end[$j]\n";
      my $window = $beg[$i]."_$end[$j]";
      my $len = $end[$j] - $beg[$i] + 1;
      while ($tempmatches[1] < $end[$j]) {
	shift @tempmatches;
	shift @tempmatches;
      }
      next if !$tempmatches[0];
      for (my $k = 0; $k <= $#tempmatches; ++$k) {
	$score{$window} += $len if $tempmatches[$k] <= $beg[$i];
      }
    }
  }
}

my $n = 0;
if ($copiesonly) {
 print "Minimum window size: $minimumwindow bp 
 Begin end length copynr \n";
} else {
  print "Minimum window size: $minimumwindow bp
Begin end length copynr score \n";
}

my @sortedbyscore = sort byscore keys %score; 
my $start = $#sortedbyscore - 50;
$start = 0 if $start < 0;
for ($i = $start; $i <=  $#sortedbyscore; ++$i) {
  my $begend =  $sortedbyscore[$i];
  my $beg = $begend;
  my $end = $beg;
  $beg =~ s/_\d+$//;
  $end =~ s/^\d+_//;
  my $len = $end - $beg + 1;
  if ($copiesonly) {
    print "$beg $end $len $score{$begend}\n";
  } else {
    my $nr = $score{$begend} / $len;
    print "$beg $end $len $nr $score{$begend}\n";
  }
  ++$n;
  last if $n == 50;
}

sub byscore {
  $score{$a} <=> $score{$b};
} 
