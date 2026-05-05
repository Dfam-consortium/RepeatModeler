#!/usr/local/bin/perl
use strict;
use warnings;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use Data::Dumper;
use File::Path qw(make_path);
use File::Copy qw(copy);
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempdir);
use File::Basename qw(basename);
use File::Spec;

use RepModelConfig;
my $RMSK_DIR = $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
my $RMBLAST_DIR = $RepModelConfig::configuration->{'RMBLAST_DIR'}->{'value'};
my $RMBLASTN = "$RMBLAST_DIR/rmblastn";
my $BLASTX = "$RMBLAST_DIR/blastx";
my $dbFile          = "$RMSK_DIR/Libraries/RepeatMasker.lib";
my $dbProtFile      = "$RMSK_DIR/Libraries/RepeatPeps.lib";
my $STKTOSAM        = "$FindBin::RealBin/stkToSam.py";
$ENV{"BLASTMAT"}    = "$FindBin::RealBin/../Matrices/ncbi/nt";

# TODO: This should be configured through RepeatModeler's configure script and located
#       in the RepModelConfig module.
my $ULTRA           = "/usr/local/ultra/ultra";
# TODO: This should be handled differently...this dependency should be cached in this project
# ......or made available through a standard URL
my $igv_source = "https://www.repeatmasker.org/~rhubley/igvnew/dist/igv.esm.min.js";

# See if there is a reasonable default for the output directory
my $defaultOutDir;
my $urlMapping;
if ( $ENV{HOME} && -d $ENV{HOME} ) {
    if ( -d $ENV{HOME} . "/public_html" && -w $ENV{HOME} . "/public_html" ) {
        $defaultOutDir = $ENV{HOME} . "/public_html/te_browser";
        ## TODO: This is specific to our server..make this more general
        $urlMapping = "https://www.repeatmasker.org/~" . basename($ENV{HOME}) . "/te_browser/index.html";
    }
}

my $public_html     = "/home/rhubley/public_html/default";
my $defaultHTML     = "$public_html/default.html";
my $defaultRef      = "$public_html/ref.fa";
my $SAMFile         = "$public_html/seed.sam";
my $MAX_TRACK_DEPTH = 10;
my @distinctColors = generate_color_palette();

### Options ###
my $outdir;
my $DEBUG = 0;

GetOptions(
    'outdir=s' => \$outdir,
    'debug!'   => \$DEBUG,
) or die "Usage: $0 [-outdir <path>] [-debug] <sequence file>|DF#########\n";

### Decide final destination ###
my $finalOutDir;
if (defined $outdir && $outdir ne '') {
    $finalOutDir = $outdir;
} elsif (defined $defaultOutDir && $defaultOutDir ne '') {
    $finalOutDir = $defaultOutDir;
} else {
    die "No output directory specified. Use -outdir <path>, or define a writable default (e.g. \$HOME/public_html/...)\n";
}

# Create a temp work directory (everything goes here first)
my $workdir = tempdir("te_browser_XXXXXX", TMPDIR => 1, CLEANUP => ($DEBUG ? 0 : 1));

# Workdir output files (published later)
my $paths = {
    workdir => $workdir,
    html    => File::Spec->catfile($workdir, "index.html"),
    ref     => File::Spec->catfile($workdir, "ref.fa"),
    sam     => File::Spec->catfile($workdir, "seed.sam"),
};


### Main Execution ###
my $seqFile = shift @ARGV or die "Usage: $0 [-outdir <path>] [-debug] <sequence file>|DF#########\n";

print "#\n";
print "# te_browser.pl\n";
print "#\n";
print "#   RMBLASTN:   $RMBLASTN\n";
print "#   BLASTX:     $BLASTX\n";
print "#   ULTRA:      $ULTRA\n";
print "#   DB (DNA):   $dbFile\n";
print "#   DB (Prot):  $dbProtFile\n";
print "#   STK->SAM:   $STKTOSAM\n";
print "#   Work dir:   $workdir\n";
print "#   Publish to: $finalOutDir\n";
print "#   Max depth:  $MAX_TRACK_DEPTH\n";
print "#   Debug:      " . ($DEBUG ? "ON" : "OFF") . "\n";
print "\n";

my ($seqID, $hasSeed, $finalSeqFile) = prepare_sequence($seqFile, $paths, $DEBUG);

print "SequenceID: $seqID\n";

# Run steps with simple counts
my $self_raw    = run_self_alignment($finalSeqFile, $DEBUG);
my $self_clust  = cluster_self_alignments($self_raw);
print "Self alignments: " . scalar(@$self_raw) . " raw, " . scalar(@$self_clust) . " after clustering/depth\n";

my $repeat_ann  = run_repeat_alignment($finalSeqFile, $DEBUG);
print "Repeat (DNA) alignments: " . scalar(@$repeat_ann) . "\n";

my $prot_ann    = run_protein_alignment($finalSeqFile, $DEBUG);
print "Protein alignments: " . scalar(@$prot_ann) . "\n";

my $ultra_ann   = run_ultra_annotation($finalSeqFile, $DEBUG);
print "ULTRA annotations: " . scalar(@$ultra_ann) . "\n\n";

my %annotations = (
    self    => $self_clust,
    repeat  => $repeat_ann,
    protein => $prot_ann,
    ultra   => $ultra_ann,
);

generate_html($seqID, $finalSeqFile, \%annotations, $hasSeed, $paths->{html}, $paths->{ref}, $paths->{sam}, $igv_source);

publish_results($paths->{workdir}, $finalOutDir);

print "Annotation complete.\n";
print "Output: $finalOutDir/index.html\n";
if ( $defaultOutDir && $urlMapping ) {
    print "URL: $urlMapping\n";
}

### Subroutines ###

sub publish_results {
    my ($workdir, $finalOutDir) = @_;
    make_path($finalOutDir) unless -d $finalOutDir;

    opendir(my $dh, $workdir) or die "Can't open workdir $workdir: $!\n";
    while (my $e = readdir($dh)) {
        next if $e eq '.' || $e eq '..';
        my $src = File::Spec->catfile($workdir, $e);
        my $dst = File::Spec->catfile($finalOutDir, $e);

        next unless -f $src;  # this script only writes files

        # Try fast move (rename), fallback to copy+unlink (cross-filesystem safe)
        if (!rename($src, $dst)) {
            copy($src, $dst) or die "Failed to copy $src -> $dst: $!\n";
            unlink($src);
        }
    }
    closedir($dh);
}

sub compute_cigar {
    my ($ref, $cons) = @_;
    my $cigar = '';
    my $count = 0;
    my $op    = '';
    my $len   = length($ref);

    for (my $i = 0; $i < $len; $i++) {
        my $r = substr($ref,  $i, 1);
        my $c = substr($cons, $i, 1);
        my $this_op;

        if ($r eq '-' && $c ne '-') {
            $this_op = 'D';  # deletion
        } elsif ($r ne '-' && $c eq '-') {
            $this_op = 'I';  # insertion 
        } else {
            $this_op = 'M';  # match/mismatch
        }

        if ($this_op eq $op) {
            $count++;
        } else {
            $cigar .= $count . $op if $op;
            $op = $this_op;
            $count = 1;
        }
    }
    $cigar .= $count . $op if $op;
    return $cigar;
}


sub cluster_self_alignments {
    my ($alignments) = @_;
    # Sort by query_start ascending, query_end descending (i.e. lowest start, longest range first)
    my @aligns = sort {
        $a->{ref_start} <=> $b->{ref_start}
          ||
        $b->{ref_end} <=> $a->{ref_end}
    } @$alignments;

    my %used_aligns;
    my @results;
    for (my $i = 0; $i < @aligns; $i++) {
        my $align = $aligns[$i];
        # Skip perfect diagonal
        next if ($align->{ref_start} == $align->{cons_start} && $align->{ref_end} == $align->{cons_end});

        my $key = join(":", $align->{ref_start}, $align->{ref_end}, $align->{cons_start}, $align->{cons_end});
        my $key_sym = join(":", $align->{cons_start}, $align->{cons_end}, $align->{ref_start}, $align->{ref_end});
        next if $used_aligns{$key} || $used_aligns{$key_sym};

        my $range_max = $align->{ref_end};
        my @cluster = ($align);
        my @ends = ($align->{ref_end});

        # Find other alignments that are close in query_start and overlap in subject
        for (my $j = 0; $j < @aligns; $j++) {
            next if $j == $i;
            my $other = $aligns[$j];

            # Only consider unused alignments
            my $other_key = join(":", $other->{ref_start}, $other->{ref_end}, $other->{cons_start}, $other->{cons_end});
            my $other_key_sym = join(":", $other->{cons_start}, $other->{cons_end}, $other->{ref_start}, $other->{ref_end});
            next if $used_aligns{$other_key} || $used_aligns{$other_key_sym};
            next if ($other->{ref_start} == $other->{cons_start} && $other->{ref_end} == $other->{cons_end});

            # Clustering condition (allowing small slop of 5 positions)
            if (
                $other->{ref_start} > $align->{ref_start} - 5 &&
                $other->{ref_start} < $align->{ref_start} + 5 &&
                $other->{cons_start} < $align->{ref_end}
            ) {
                if ($other->{ref_end} > $range_max) {
                    $range_max = $other->{ref_end};
                }
                push @cluster, $other;
                push @ends, $other->{ref_end};
            }
        }

        if (@cluster > 2) {
            # Collapse cluster into a single alignment
            my $max_query = $align->{ref_end};
            foreach my $other (@cluster) {
                my $okey = join(":", $other->{ref_start}, $other->{ref_end}, $other->{cons_start}, $other->{cons_end});
                my $okey_sym = join(":", $other->{cons_start}, $other->{cons_end}, $other->{ref_start}, $other->{ref_end});
                $used_aligns{$okey} = 1;
                $used_aligns{$okey_sym} = 1;
                $max_query = $other->{ref_end} if $other->{ref_end} > $max_query;
            }
            push @results, {
                ref_start    => $align->{ref_start},
                ref_end      => $max_query,
                cons_start   => $align->{ref_start},
                cons_end     => $max_query,
                orient       => "+",
                name         => "tandem_cluster",
                score        => $align->{score},
                type         => "tandem",
                ends         => [@ends]
            };
        } else {
            foreach my $other (@cluster) {
                my $okey = join(":", $other->{ref_start}, $other->{ref_end}, $other->{cons_start}, $other->{cons_end});
                my $okey_sym = join(":", $other->{cons_start}, $other->{cons_end}, $other->{ref_start}, $other->{ref_end});
                $used_aligns{$okey} = 1;
                $used_aligns{$okey_sym} = 1;
                push @results, $other;
            }
        }
    }
    return \@results;
}

sub read_sequence_ids {
    my ($file) = @_;
    open my $in, "fgrep '>' $file |" or die "Can't open $file\n";
    my @ids = map { /^>(\S+)/ ? $1 : () } <$in>;
    close $in;
    return @ids;
}

sub run_self_alignment {
    my ($file, $DEBUG) = @_;
    my $cmd = "$RMBLASTN " . get_blast_params('self') . " -subject $file -query $file";
    print "CMD: $cmd\n" if $DEBUG;
    return parse_blast_output($cmd, 'self', $DEBUG);
}

sub run_repeat_alignment {
    my ($file, $DEBUG) = @_;
    my $cmd = "$RMBLASTN " . get_blast_params('repeat') . " -db $dbFile -query $file";
    print "CMD: $cmd\n" if $DEBUG;
    return parse_blast_output($cmd, 'repeat', $DEBUG);
}

sub run_protein_alignment {
    my ($file, $DEBUG) = @_;
    my $cmd = "$BLASTX -word_size 2 -evalue 0.001 -db $dbProtFile -query $file -outfmt \"6 evalue perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand sseqid sstart send slen qframe sseq\"";
    print "CMD: $cmd\n" if $DEBUG;
    return parse_blast_output($cmd, 'protein', $DEBUG);
}

sub run_ultra_annotation {
    my ($file, $DEBUG) = @_;
    my $cmd = "$ULTRA $file";
    print "CMD: $cmd\n" if $DEBUG;
    open my $in, "$cmd |" or die "Can't run $cmd\n";
    my @annotations;
    while (<$in>) {
        next unless /^\S+\t\d+\t/;
        my @data = split(/\t/);
        my $name = $data[3];
        # If name is longer than 10 bp, shorten as 4 + '..' + 4
        if (length($name) > 14) {
            $name = substr($name, 0, 4) . ".." . substr($name, -4) . "[" . length($name) . "]";
        }
        push @annotations, {
            ref_start => $data[1],
            ref_end   => $data[2],
            name      => $name,
            color     => "blue",
            type      => 'ultra'
        };
    }
    close $in;
    return \@annotations;
}



sub parse_blast_output {
    my ($cmd, $type, $DEBUG) = @_;
    open my $in, "$cmd |" or die "Can't run $cmd\n";
    my @raw_annots;

    while (<$in>) {
        chomp;
        my @d = split(/\t/);

        my ($ref_start, $ref_end, $cons_start, $cons_end, $orient, $name, $score, $ref_seq, $cons_seq, $oseq, $cons_len);

        if ($type eq 'protein') {
            # BLASTX output:
            # 0: evalue
            # 1: perc_sub
            # 2: perc_query_gap
            # 3: perc_db_gap
            # 4: qseqid
            # 5: qstart
            # 6: qend
            # 7: qlen
            # 8: qstrand
            # 9: sseqid
            # 10: sstart
            # 11: send
            # 12: slen
            # 13: qframe
            # 14: sseq
            $score      = $d[0];
            if ( $d[5] <= $d[6] ) { 
              $orient = "+";
              $ref_start  = $d[5];
              $ref_end    = $d[6];
            }else {
              $orient = "-";
              $ref_start  = $d[6];
              $ref_end    = $d[5];
            }
            $cons_start = $d[10];
            $cons_end   = $d[11];
            $cons_len   = $d[12];
            #$orient     = ($d[8] && $d[8] eq 'minus') ? '-' : '+';
            $name       = $d[9] || 'protein';
            $ref_seq    = undef;           # BLASTX does NOT report aligned qseq as a field
            $cons_seq   = $d[14];
            $oseq       = $d[14];
        } else {
            # RMBLASTN output:
            # 0: score
            # 1: perc_sub
            # 2: perc_query_gap
            # 3: perc_db_gap
            # 4: qseqid
            # 5: qstart
            # 6: qend
            # 7: qlen
            # 8: sstrand
            # 9: sseqid
            # 10: sstart
            # 11: send
            # 12: slen
            # 13: kdiv
            # 14: cpg_kdiv
            # 15: transi
            # 16: transv
            # 17: cpg_sites
            # 18: qseq
            # 19: sseq
            $score      = $d[0];
            $ref_start  = $d[5];
            $ref_end    = $d[6];
            $cons_start = $d[10];
            $cons_end   = $d[11];
            $cons_len   = $d[12];
            $orient     = ($d[8] && $d[8] eq 'minus') ? '-' : '+';
            $name       = $d[9] || 'self';
            $ref_seq    = $d[18];
            $cons_seq   = $d[19];
            $oseq       = $cons_seq;
        }

        # skip perfect diagonal for self
        next if ($type eq 'self' && $ref_start == $cons_start && $ref_end == $cons_end);

        push @raw_annots, {
            ref_start  => $ref_start,
            ref_end    => $ref_end,
            cons_start => $cons_start,
            cons_end   => $cons_end,
            cons_len   => $cons_len,
            orient     => $orient,
            name       => $name,
            score      => $score,
            ref_seq    => $ref_seq,
            cons_seq   => $cons_seq,
            oseq       => $oseq,
            type       => $type
        };
    }
    close $in;
    if ( $DEBUG ) {
      print "" . Dumper(\@raw_annots) . "\n";
    }
    #  print "Total Annotations: " . scalar(@raw_annots) . "\n";
    return apply_depth_limit(\@raw_annots, $type);
}


sub apply_depth_limit {
    my ($annotations, $type) = @_;
    my $max_depth = $MAX_TRACK_DEPTH;

    my @sorted = sort {
        ($type eq 'protein') ? $a->{score} <=> $b->{score} : $b->{score} <=> $a->{score}
    } @$annotations;

    my @accepted;
    my @depth_bins;

    foreach my $ann (@sorted) {
        my ($start, $end) = @{$ann}{qw(ref_start ref_end)};
        my $max_current_depth = 0;

        for my $pos ($start .. $end) {
            $max_current_depth = $depth_bins[$pos] // 0 if defined $depth_bins[$pos] && $depth_bins[$pos] > $max_current_depth;
        }

        if ($max_current_depth < $max_depth) {
            push @accepted, $ann;
            $depth_bins[$_]++ for $start .. $end;
        }
    }
    return \@accepted;
}

sub get_blast_params {
    my ($type) = @_;
    my $mask_level = ($type eq 'self') ? 101 : 80;

    return "-num_alignments 9999999 -gapopen 20 -gapextend 5 " .
           "-mask_level $mask_level -complexity_adjust -word_size 14 " .
           "-xdrop_ungap 400 -xdrop_gap_final 200 -xdrop_gap 100 " .
           "-min_raw_gapped_score 200 -dust no " .
           "-outfmt \"6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand " .
           "sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites qseq sseq\" -matrix comparison.matrix";
}

sub chain_alignments {
    my ($annotations, $max_gap, $max_overlap) = @_;
    $max_gap     //= 100;
    $max_overlap //= 20;

    my @sorted = sort {
        $a->{name} cmp $b->{name} ||
        $a->{orient} cmp $b->{orient} ||
        $a->{ref_start} <=> $b->{ref_start}
    } @$annotations;

    my @chains;
    my $current_chain;

    foreach my $ann (@sorted) {
        if (!$current_chain 
            || $current_chain->{name} ne $ann->{name} 
            || $current_chain->{orient} ne $ann->{orient}
        ) {
            push @chains, finalize_chain($current_chain) if $current_chain;
            $current_chain = start_new_chain($ann);
            next;
        }

        # Verify colinearity in both ref and cons space
        my $ref_gap  = $ann->{ref_start} - $current_chain->{last_ref_end};
        my $cons_gap = $ann->{cons_start} - $current_chain->{last_cons_end};

        if (
            ($ref_gap >= -$max_overlap && $ref_gap <= $max_gap) &&
            ($cons_gap >= -$max_overlap && $cons_gap <= $max_gap) &&
            (($ann->{ref_start} >= $current_chain->{last_ref_end} && $ann->{cons_start} >= $current_chain->{last_cons_end}) ||
             ($ann->{ref_start} <= $current_chain->{last_ref_end} && $ann->{cons_start} <= $current_chain->{last_cons_end}))
        ) {
            # Colinear, append to current chain
            $current_chain->{ref_end} = $ann->{ref_end};
            $current_chain->{last_ref_end}  = $ann->{ref_end};
            $current_chain->{last_cons_end} = $ann->{cons_end};

            my $oseq   = $ann->{ref_seq}  // '';
            my $seq  = $ann->{cons_seq} // '';
            my $cigar = ($seq ne '' && $oseq ne '') ? compute_cigar($seq, $oseq) : '';
            push @{ $current_chain->{components} }, { 
                start  => $ann->{ref_start}, 
                end    => $ann->{ref_end}, 
                ostart => $ann->{cons_start}, 
                oend   => $ann->{cons_end},
                osize  => $ann->{cons_len},
                seq    => ($seq  =~ tr/-//dr),
                oseq   => ($oseq =~ tr/-//dr),
                cigar  => $cigar
            };

        } else {
            # Start new chain
            push @chains, finalize_chain($current_chain);
            $current_chain = start_new_chain($ann);
        }
    }
    push @chains, finalize_chain($current_chain) if $current_chain;
    return \@chains;
}

sub start_new_chain {
    my ($ann) = @_;
    my $oseq   = $ann->{ref_seq}  // '';
    my $seq  = $ann->{cons_seq} // '';
    my $cigar = ($seq ne '' && $oseq ne '') ? compute_cigar($seq, $oseq) : '';

    return {
        name        => $ann->{name},
        orient      => $ann->{orient},
        color       => $ann->{color} || 'gray',
        ref_start   => $ann->{ref_start},
        ref_end     => $ann->{ref_end},
        osize       => $ann->{cons_len},
        last_ref_end  => $ann->{ref_end},
        last_cons_end => $ann->{cons_end},
        components  => [ { 
            start  => $ann->{ref_start}, 
            end    => $ann->{ref_end}, 
            ostart => $ann->{cons_start}, 
            oend   => $ann->{cons_end},
            seq    => ($seq  =~ tr/-//dr),   # gapless
            oseq   => ($oseq =~ tr/-//dr),   # gapless
            cigar  => $cigar
        } ]
    };
}

sub finalize_chain {
    my ($chain) = @_;
    return unless $chain; # Fix for undefined chain
    my @starts  = map { $_->{ostart} } @{ $chain->{components} };
    my @ends    = map { $_->{oend} } @{ $chain->{components} };

    my $min_ostart = (sort { $a <=> $b } @starts)[0];
    my $max_oend   = (sort { $b <=> $a } @ends)[0];

    $chain->{ostart} = $min_ostart;
    $chain->{oend}   = $max_oend;

    return $chain;
}

# Helper to convert 1-based, inclusive to 0-based, half-open
sub to_0_based_half_open {
    my ($start, $end) = @_;
    return ($start - 1, $end); # $end is already inclusive, so leave as-is for output, but must add 1 for half-open in JS
}




sub prepare_sequence {
    my ($file, $paths, $DEBUG) = @_;
    my $hasSeed = 0;
    my $finalSeqFile = $file;

    my $stkFile;
    if ( ! -s $file && $file =~ /^D[FR]\d{9}$/ ) {
        print "Source: Dfam family ($file), downloading seed alignment from:\n";
        print "     https://www.dfam.org/api/families/$file/seed?format=stockholm\n";
        system("wget \"https://www.dfam.org/api/families/$file/seed?format=stockholm\" -O tmpAnnotSeqDfamSeed.stk > /dev/null 2>&1");
        die "Failed to download seed for $file\n" unless -s "tmpAnnotSeqDfamSeed.stk";

       $stkFile = "tmpAnnotSeqDfamSeed.stk";
    }elsif ( $file =~ /\.stk$/ ) {
        print "Source: Stockholm file ($file)\n";
        $stkFile = $file;
    }

    if ( $stkFile ) {
        open STK,"<$stkFile" or die "Can't open $stkFile\n";
        my $RF = "";
        while (<STK>) {
          if ( /^#=GC\s+RF\s+(\S+)/ ) {
            $RF = $1;
            last;
          }
        }
        close STK;
        if ( $RF =~ /^[\.Xx]+$/ ) {
          print "  ** Stockholm file uses occupancy-only RF line (e.g 'x's).  Rebuilding consensus sequence...\n";
          # No consensus, must run Linup
          system("/home/rhubley/projects/RepeatModeler/util/Linup -stockholm $stkFile > tmpFixed.stk");
          system("mv tmpFixed.stk tmpAnnotSeqDfamSeed.stk");
          $stkFile = "tmpAnnotSeqDfamSeed.stk";
        }
 
        print "Converting to SAM format for IGV visualization...\n";
        system("$STKTOSAM $stkFile > stkToSam.log");

        if (-s "out.sam" && -s "ref.fa") {
            rename("ref.fa", "tmpAnnotSeqDfamCons.fa");
            system("mv out.sam $paths->{sam}");
            $finalSeqFile = "tmpAnnotSeqDfamCons.fa";
            $hasSeed = 1;
        } else {
            die "Failed to process seed alignment for accession $file.\n";
        }
    }

    my @ids = read_sequence_ids($finalSeqFile);
    die "Seqfile $finalSeqFile contains more than one sequence!\n" if @ids > 1;

    copy($finalSeqFile, $paths->{ref});
    return ($ids[0], $hasSeed, $finalSeqFile);
}

sub generate_html {
    my ($seqID, $seqFile, $annotsRef, $hasSeed, $outFile, $refFile, $samFile, $igv_js_url) = @_;
    print "TempDir: $outFile\n";
    open my $out, "> $outFile" or die "Can't write to $outFile\n";

    print $out <<"HTML_HEADER";
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>$seqID Annotation</title>
</head>
<body>
    <h1>$seqID Annotation</h1>
    <div id="igvDiv" style="padding:10px; border:1px solid lightgray;"></div>

    <script type="module">
        import igv from "$igv_js_url";

        const options = {
            reference: { fastaURL: "ref.fa", indexed: false },
            locus: ["$seqID"],
            tracks: [
HTML_HEADER

    if ($hasSeed) {
        print $out <<"SEED_TRACK";
                {
                    name: 'Seed Alignment',
                    type: 'seedalign',
                    format: 'sam',
                    sourceType: 'sam',
                    url: "seed.sam",
                    displayMode: "SQUISHED",
                    autoHeight: true
                },

SEED_TRACK
    }

    for my $type (qw/ultra self/) {
        emit_track($out, $seqID, $type, $annotsRef->{$type});
    }

    my $repeat_chains = chain_alignments($annotsRef->{repeat});
    emit_chain_track($out, $seqID, 'TE DNA Homology', $repeat_chains);

    emit_track($out, $seqID, 'protein', $annotsRef->{'protein'});
    print $out <<"HTML_FOOTER";
            ]
        };

        igv.createBrowser(document.getElementById("igvDiv"), options)
            .then(browser => console.log("IGV browser created."))
            .catch(error => console.error("Error creating IGV browser:", error));
    </script>
</body>
</html>
HTML_FOOTER

    close $out;
}
sub emit_chain_track {
    my ($out, $seqID, $trackName, $chains) = @_;
    return unless @$chains;

    print $out <<"CHAIN_START";
                {
                    name: '$trackName',
                    type: 'chain',
                    displayMode: 'EXPANDED',
                    autoHeight: true,
                    features: [
CHAIN_START

    my %color_by_subject = ();
    my $color_idx = 0;
    foreach my $chain (@$chains) {
        my $components_json = join(",\n", map {
            "      { start: ".($_->{start}-1).", end: ".($_->{end}+0).", ostart: ".($_->{ostart}-1).", oend: ".($_->{oend}+0).", seq: '".($_->{seq}//'')."', oseq: '".($_->{oseq}//'')."', cigar: '".($_->{cigar}//'')."' }"
        } @{ $chain->{components} });

        my $subject_id = $chain->{name};
        if (!exists $color_by_subject{$subject_id}) {
            $color_by_subject{$subject_id} = $distinctColors[$color_idx++ % @distinctColors];
        }
        my $color = $color_by_subject{$subject_id};

        print $out <<"FEATURE";
                        {
                            chr: "$seqID",
                            start: @{[ $chain->{ref_start} - 1 ]},
                            end: @{[ $chain->{ref_end} + 0 ]},
                            name: "$chain->{name}",
                            strand: "$chain->{orient}",
                            color: "$color",
                            ostart: @{[ $chain->{ostart} - 1 ]},
                            oend: @{[ $chain->{oend} + 0 ]},
                            osize: $chain->{osize},
                            components: [
$components_json
                            ]
                        },
FEATURE
    }

    print $out <<"CHAIN_END";
                    ]
                },
CHAIN_END
}

sub emit_track {
    my ($out, $seqID, $trackName, $annotations) = @_;
    return unless @$annotations;

    my $trackLabel = {
        ultra   => 'ULTRA Annotations',
        self    => 'Self Alignments',
        repeat  => 'TE DNA Homology',
        protein => 'TE Protein Homology'
    }->{$trackName} || ucfirst($trackName);

    # Special handling for "protein" track in UCSC PSL-like style
    if ($trackName eq 'protein') {
        print $out <<"TRACK_START";
                {
                    name: "TE Protein Homology",
                    type: 'annotation',
                    autoHeight: true,
                    features: [
TRACK_START

        my %color_by_subject;
        my $color_idx = 0;
        foreach my $annot (@$annotations) {
            my $subject_id = $annot->{name};
            my $color = exists $color_by_subject{$subject_id}
                ? $color_by_subject{$subject_id}
                : ($color_by_subject{$subject_id} = $distinctColors[$color_idx++ % @distinctColors]);

            # Normalize coordinates and strand
            my $strand = (defined $annot->{orient} && $annot->{orient} eq '-') ? '-' : '+';

            my ($start, $end) = ($annot->{ref_start} - 1, $annot->{ref_end} + 0);
            if ($start > $end) {
                ($start, $end) = ($end, $start);
            }
            my ($cdStart, $cdEnd) = ($start, $end);

            my $exon_json = "{start:$start, end:$end, cdStart:$cdStart, cdEnd:$cdEnd}";

            my $oChromStart = ($annot->{cons_start} || 1) - 1;
            my $oChromEnd   = $annot->{cons_end} || ($oChromStart+1);
            my $oStrand     = (defined $annot->{orient} && $annot->{orient} ne '') ? $annot->{orient} : "+";
            my $oChromSize  = abs($oChromEnd - $oChromStart);

            my $oChromStarts = "$oChromStart,";
            my $oSequence = $annot->{oseq} // '';

            print $out <<"PROTEIN_FEATURE";
                            {
                                "chr":"$seqID",
                                "start":$start,
                                "end":$end,
                                "name":"$subject_id",
                                "score":@{[$annot->{score}||0]},
                                "strand":"$strand",
                                "cdStart":$cdStart,
                                "cdEnd":$cdEnd,
                                "color":"$color",
                                "exons":[$exon_json],
                                "oChromStart":"$oChromStart",
                                "oChromEnd":"$oChromEnd",
                                "oStrand":"$oStrand",
                                "oChromSize":"$oChromSize",
                                "oChromStarts":"$oChromStarts",
                                "oSequence":"$oSequence"
                            },
PROTEIN_FEATURE
        }
        print $out <<"TRACK_END";
                    ]
                },
TRACK_END
        return;
    }

    if ( $trackName eq 'self' ) {
        # Special case for self alignments
        print $out <<"SELF_TRACK_START";
                {
                    name: '$trackLabel',
                    type: 'selfpair',
                    autoHeight: true,
                    features: [
SELF_TRACK_START
    }else {
        print $out <<"OTHER_TRACK_START";
                {
                    name: '$trackLabel',
                    type: 'annotation',
                    autoHeight: true,
                    format: 'bed',
                    features: [
OTHER_TRACK_START
    }

    my %color_by_subject;
    my $color_idx = 0;

    foreach my $annot (@$annotations) {
        my $subject_id = $annot->{name};
        my $color;
        if ($trackName eq 'repeat' || $trackName eq 'protein') {
            if (!exists $color_by_subject{$subject_id}) {
                $color_by_subject{$subject_id} = $distinctColors[$color_idx++ % @distinctColors];
            }
            $color = $color_by_subject{$subject_id};
        } else {
            $color = $annot->{color} || $distinctColors[$color_idx++ % @distinctColors];
        }

        # SPECIAL CASE: self alignments
        if ($trackName eq 'self') {
               my $minStart = $annot->{ref_start} - 1;
               $minStart = $annot->{cons_start} - 1 if $annot->{cons_start} < $minStart;
               my $maxEnd = $annot->{ref_end} + 0;
               $maxEnd = $annot->{cons_end} + 0 if $annot->{cons_end} > $maxEnd;

               my $minCons = $annot->{cons_start} - 1;
               my $maxCons = $annot->{cons_end} + 0;
               if ( $annot->{cons_start} - 1 > $annot->{cons_end} + 0 ) {
                 $minCons = $annot->{cons_end} + 0;
                 $maxCons = $annot->{cons_start} - 1;
               }
               print $out <<"SELF_PAIR";
                {
                    chr: "$seqID",
                    start: $minStart,
                    end: $maxEnd,
                    pstart: @{[ $annot->{ref_start} - 1 ]},
                    pend: @{[ $annot->{ref_end} + 0 ]},
                    sstart: $minCons,
                    send: $maxCons,
                    color: "$color",
                    strand: "$annot->{orient}" 
                },
SELF_PAIR
               #}
        } else {
            # Non-self: one feature as before
            print $out <<"FEATURE_START";
                        {
                            chr: "$seqID",
                            start: @{[ $annot->{ref_start} - 1 ]},
                            end: @{[ $annot->{ref_end} + 0 ]},
                            color: "$color",
FEATURE_START

            if (defined $annot->{orient} && $annot->{orient} ne '') {
                print $out "                            strand: \"$annot->{orient}\",\n";
            }

            print $out <<"FEATURE_END";
                            name: "$annot->{name}"
                        },
FEATURE_END
        }
    }

    print $out <<"TRACK_END";
                    ]
                },
TRACK_END
}

sub generate_color_palette {
return ('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000', '#ff9180', '#a65b29', '#665a4d', '#998226', '#143300', '#7fff80', '#6cd9d2', '#0077b3', '#265499', '#1f00e6', '#6c468c', '#d9a3d5', '#a6296c', '#d96c7b', '#4c0000', '#8c5946', '#664733', '#bf8000', '#736739', '#74a653', '#004d1f', '#00474d', '#99bbcc', '#80b3ff', '#0a004d', '#b63df2', '#f200c2', '#997387', '#66333a', '#330d0d', '#ffd0bf', '#ff8800', '#332b1a', '#d9d26c', '#299900', '#008044', '#46888c', '#102940', '#293aa6', '#290099', '#2b1a33', '#e673cf', '#660029', '#d9a3aa', '#806060', '#332a26', '#331b00', '#f2deb6', '#add900', '#d0ffbf', '#73e6b0', '#bffbff', '#406280', '#7382e6', '#2e1966', '#8f0099', '#594355', '#ff80b3', '#e55039', '#ff6600', '#ffc480', '#f2c200', '#475900', '#53664d', '#0d3326', '#0099bf', '#0061f2', '#7c82a6', '#c6b6f2', '#300033', '#f20081', '#f20041', '#66241a', '#592400', '#a67f53', '#4c3d00', '#aab386', '#20f200', '#468c75', '#39494d', '#001f4d', '#4d5066', '#a173e6', '#591655', '#400022', '#990014' );

}
