#!/usr/local/bin/perl 
##---------------------------------------------------------------------------##
##  File:
##      @(#) align.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Arian Smit         asmit@systemsbiology.org
##  Description:
##      A wrapper script around cross_match, rmblastn and
##      nhmmer with a common command line interface.
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
=head1 NAME

align.pl - A standard interface for several sequence aligners

=head1 SYNOPSIS

  align.pl [-version]
or
  align.pl [options] <fasta seq1> <fasta seq2>
or
  align.pl [options] -query <fasta seq1> -database (<fasta seq2>)

  TBD: align.pl [options] <fasta seq1> (<fasta seq2>)

  options            
  -alignments,-a    Instruct aligner to produce full alignment details. (default: off)
  -bandwidth,-ba    For aligners that use a banded smith-waterman approach, this parameter
                      sets the bandwidth.  This value is interpreted to be the max offset
                      from the diagonal, the full bandwidth is 2*bandwidth+1. (default: 14)
  -gap_init,-g      The gap open penalty for either an insertion or deletion. (default: -25)
                      The convention used by this interface (and cross_match) is that the
                      gap_init penalty is equivalent to a gap initiation of one residue all
                      subsequent gap residues are scored using the extension penalties.
  -extension,-e     The general gap extension penalty for aligners that do not support
                      different penalties for insertions/deletions.  This is the equivalent
                      of the cross_match -gap_ext parameter. (default: -5)
  -del_gap_ext,-d   The penalty for a deletion gap extension for aligners that support
                      different penalties for insertions/deletions.  If this option is set
                      without -ins_gap_ext the -ins_gap_ext penalty is automatically set to 
                      (del_gap_ext - 1). For aligners that do not support this distinction,
                      del_gap_ext is used for all extensions.  
  -ins_gap_ext,-i   The penalty for an insertion gap extension for aligners that support
                      different penalties for insertions/deletions. See del_gap_ext.
  -masklevel,-l     Filter alignments by query overlap.  -masklevel (default: 80)
  -matrix,-ma       Scoring matrix for residue substitution.  Resolution of a path to the
                      matrix is handled by first checking if the matrix path resolves to a
                      file directly.  If not, the environment variable MATRIX_DIR is 
                      checked for the existence of the named matrix.  Finally, the RepeatModeler
                      and RepeatMasker configured matrix directories are checked in that order.  
                      Aliases for common RepeatMasker/RepeatModeler matrices are further supported:
                         - A single number resolves to "##p41g.matrix" 
                         - A string "##p##" resolves to "##p##g.matrix"
  -minmatch,-mi     The minimum length alignment to report. (default: 14)
  -raw,-r           -raw (default off, i.e. score is complexity adjusted) 
  -minscore,-s      -minscore (default 200)
  -word_raw,-w      -word_raw (default off, i.e. length of seeds is higher with lower complexity)
  -x                This invokes the '-screen' option if the crossmatch search engine is
                    used. This will create a *.screen file where all aligned regions of the
                    query are replaced with "X"s. (default off)
  -threads,-t       Set the number of threads to use for rmblastn (cross_match is single threaded).
                    (default: 6)
  -mt_qmode,-mtq    Set the threading scheme for rmblastn (2.13+) to query mode (mt_mode 1). The
                    default rmblastn scheme is to thread by subject. (default: off)
  -crossmatch,-cm   Use cross_match search engine.  Otherwise use the default search
                      engine (configured at installation).
  -rmblast,-rm      Use RMBlast as the search engine. Otherwise use the default search
                      engine (configured at installation).
  -nhmmer,-nh       Use nhmmer as the search engine. Otherwise use the default search
                      engine (configured at installation). 
  -fmindex,-fm      If using the nhmmer search engine this will turn on the use of the
                      fmindex acceleration.
  -force,-fo        When using RMBlast this option will force the rebuilding of the frozen
                      database files, if present.
  -caf              Produce CAF output rather than standard cross_match formatting (this
                      option implies -alignments will be used).
  -bed3             Produce BED3 format: sequence_id<tab>seq_start<tab>seq_end
  -bed6             Produce BED6 format: sequence_id<tab>seq_start<tab>seq_end<tab>q_name<tab>score<tab>strand
  -k_param,-k       ALP K parameter for E_value calculation
  -lambda,-la       ALP lambda parameter for E_value calculation
  -q_size           Required when -k_param and -lambda are provided for evalue calculation

  Xmatch.pl options not yet supported:
  #-original,-o      keep the original cross_match mismatch level (mismatches/length_of_query) 
  #                  By default align.pl recalculates the mismatch level as mismatches/aligned_bases
  #-quiet,-q         STDERR of cross_match gets redirected to a file \"tempxmatch.stderr\" 
  #
  #-nab,-n           Only prints the summary lines and adds a \"+\" for forward match to give constant nr of columns
  #-zip <xmoutput>   Skip the cross_match part and take the indicated cross_match output or Xmatch.pl .rawxmout output
  #                  (e.g. to increase the minimum score displayed (-score) or to add or skip -perbasescore, -cg, or -nab)
  #
  #   # Not reliable with -blast option, until rmblast.pl has a -raw option: 
  #-perbasescore,-p  Adds a first column reporting the raw score per base of the query 
  #                 NOTE1: this forces the -raw option, so only use this on established matches, and the -cg option
  #                 NOTE2:  unless a -matrix option is used, this mode will use ~asmit/Matrices/symmetricNscores0matrix 
  #   #  For now only works within a -p run.
  #-cg              Only works and effects a -p run right for now. Subtracts from mismatch level and adjusts up scores for
  #                 one transition/ambiguous mismatch per CG in either consensus and treats TG/CA pairs as one transition.
  #                 Forces -a (as alignments are needed for the adjustment)

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2022-2023 Robert Hubley, Institute for Systems Biology

=head1 LICENSE

This code may be used in accordance with the Open Source License v. 3.0
http://opensource.org/licenses/OSL-3.0

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>
Arian Smit <asmit@systemsbiology.org>

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
use HMMERSearchEngine;
use SearchResult;
use SearchResultCollection;
use RepeatMaskerConfig;

my $Version = "0.1";
my $DEBUG = 0;

#
# Paths
#
my $CM_DIR = $RepeatMaskerConfig::configuration->{'CROSSMATCH_DIR'}->{'value'};
my $RMSK_DIR = $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
my $RMBLAST_DIR = $RepModelConfig::configuration->{'RMBLAST_DIR'}->{'value'};
my $defaultEngine = "rmblast";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    'version|v', # print out the version and exit
    'alignments|a',
    'bandwidth|ba=s',
    'blast|bl',
    'caf',
    'bed3',
    'bed6',
    'crossmatch|cm',
    'cg',
    'database=s',
    'del_gap_ext|d=i',
    'ins_gap_ext|i=i',
    'k_param|k=s',
    'lambda|la=s',
    'q_size=s',
    'extension|e=i',
    'fmindex|fm',
    'force|fo',
    'gap_init|g=i',
    'masklevel|level|l=i',
    'matrix|ma=s',
    'minmatch|mi=i',
    'minscore|score|s=i',
    'mt_qmode|mtq',
    'nab|n',
    'nhmmer|nh',
    'original|o',
    'perbasescore|p',
    'quiet|q',
    'query=s',
    'verbose|verb',
    'threads|t=i',
    'raw|r',
    'rmblast|rm',
    'rmb_db_softmask',
    'word_raw|w',
    'screen|x',
    'zip=s'
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ($options{'version'}) {
  print "$Version\n";
  exit;
}

my $queryFile;
my $databaseFile;
if ( -e $options{'query'} && -e $options{'database'} ) {
  $queryFile = $options{'query'};
  $databaseFile = $options{'database'};
}elsif ( scalar(@ARGV) == 2 ) {
  $queryFile = $ARGV[0];
  $databaseFile = $ARGV[1];
}else {
  print "\n\nMissing either '-query <file> -database <file>' or '<query_file> <database_file>'\n\n";
  usage();
}

my $engine = $defaultEngine;
$engine = "crossmatch" if ( exists $options{'crossmatch'} );
$engine = "rmblast" if ( exists $options{'rmbalst'} );
$engine = "nhmmer" if ( exists $options{'nhmmer'} );

my $gen_alignments = 0;
if ( $options{'alignments'} || $options{'caf'} ) {
  $gen_alignments = 1;
}

my $engine_dir;
my $engine_prg;
my $sEngineObj;
if ( $engine eq "crossmatch" )
{
  $engine = "crossmatch";
  unless ( -d $CM_DIR && -x "$CM_DIR/cross_match") {
    # fall back to path resolution
    my $retVal = `whereis cross_match`;
    if ( $retVal =~ /^cross_match:\s+(\S+)/ ) {
      $engine_prg = $1;
      $engine_dir = dirname($engine_prg);
    }else {
      print "\n\nERROR: Could not locate the cross_match program.  Perhaps, the configured RepeatMasker\n" .
            "       program ($RMSK_DIR) is not configured to use that\n" .
            "       search engine.\n\n";
      exit;
    }
  }else {
    $engine_dir = $CM_DIR;
    $engine_prg = "$CM_DIR/cross_match";
  }

  $sEngineObj = CrossmatchSearchEngine->new( pathToEngine => $engine_prg );
  my $params = "";
  $params .= " -screen " if ( $options{'screen'} );
  # TODO: Add to SearchEngineI and CrossMatchSearchEngine
  #$sEngineObj->setAdditionalParameters($params);

}elsif ( $engine eq "rmblast" ) {
  $engine = "rmblast";
  unless ( -d $RMBLAST_DIR && -x "$RMBLAST_DIR/rmblastn") { 
    # fall back to path resolution
    my $retVal = `whereis rmblastn`;
    if ( $retVal =~ /^rmblastn:\s+(\S+)/ ) {
      $engine_prg = $1;
      $engine_dir = dirname($engine_prg);
    }else {
      print "\n\nERROR: Could not locate the rmblastn program.  Perhaps, the configured RMBlast\n" .
            "       program ($RMBLAST_DIR) is not installed correctly.\n\n";
      exit;
    }
  }else {
    $engine_dir = $RMBLAST_DIR;
    $engine_prg = "$RMBLAST_DIR/rmblastn"; 
  }

  $sEngineObj = NCBIBlastSearchEngine->new( pathToEngine => $engine_prg );

  # Engine specific parameters
  $sEngineObj->setCores( $options{'threads'} ? $options{'threads'} : undef );
  if ( exists $options{'threads'} && exists $options{'mt_qmode'} ) {
    $sEngineObj->setThreadByQuery(1);
  }
  if ( $options{'rmb_db_softmask'} ) {
    # See codes in section that build the database
    $sEngineObj->setAdditionalParameters("-db_soft_mask 100");
  }
}elsif ( $engine eq "nhmmer" ) {
  $engine = "nhmmer";
  my $HM_DIR = $RepeatMaskerConfig::configuration->{'HMMER_DIR'}->{'value'};
  unless ( -d $HM_DIR && -x "$HM_DIR/nhmmer") {
    # fall back to path resolution
    my $retVal = `whereis nhmmer`;
    if ( $retVal =~ /^nhmmer:\s+(\S+)/ ) {
      $engine_prg = $1;
      $engine_dir = dirname($engine_prg);
    }else {
      print "\n\nERROR: Could not locate the nhmmer program.  Perhaps, the configured RepeatMasker\n" .
            "       program ($RMSK_DIR) is not configured to use that\n" .
            "       search engine.\n\n";
      exit;
    }
  }else {
    $engine_dir = $HM_DIR;
    $engine_prg = "$HM_DIR/nhmmer";
  }

  $sEngineObj = HMMERSearchEngine->new( pathToEngine => $engine_prg );
  # Engine specific parameters
  $sEngineObj->setCores( $options{'threads'} ? $options{'threads'} : undef );

}else {
  print "\n\nERROR: Could not locate the default aligner ($engine)!\n";
  exit;
}

my $engine_ver = $sEngineObj->getVersion();

  
$sEngineObj->setQuery($queryFile);
$sEngineObj->setSubject($databaseFile);
$sEngineObj->setMinScore( $options{'minscore'} ? $options{'minscore'} : 200 );
$sEngineObj->setBandwidth( $options{'bandwidth'} ? $options{'bandwidth'} : 14 );
if ( $options{'masklevel'} ne "" ) {
  $sEngineObj->setMaskLevel( $options{'masklevel'} );
}else {
  $sEngineObj->setMaskLevel( 80 );
}
$sEngineObj->setMinMatch( $options{'minmatch'} ? $options{'minmatch'} : 14 );
$sEngineObj->setGenerateAlignments( $gen_alignments );


my $gap_open_penalty = -25;
my $ins_extn_penalty = -5;
my $del_extn_penalty = -5;
if ( exists $options{'gap_init'} ){
  $gap_open_penalty = $options{'gap_init'};
}
if ( exists $options{'extension'} ) 
{
  $ins_extn_penalty = $options{'extension'};
  $del_extn_penalty = $options{'extension'};
}elsif ( exists $options{'del_gap_ext'} ) 
{
  $del_extn_penalty = $options{'del_gap_ext'};
  if ( exists $options{'ins_gap_ext'} ) {
    $ins_extn_penalty = $options{'ins_gap_ext'};
  }else {
    $ins_extn_penalty = -(abs($del_extn_penalty)+1);
  }
}elsif ( exists $options{'ins_gap_ext'} ) {
  $ins_extn_penalty = $options{'ins_gap_ext'};
  if ( exists $options{'del_gap_ext'} ) {
    $del_extn_penalty = $options{'del_gap_ext'};
  }else {
    $del_extn_penalty = -(abs($ins_extn_penalty)-1);
  }
}
$sEngineObj->setGapInit( $gap_open_penalty );
$sEngineObj->setInsGapExt( $ins_extn_penalty );
$sEngineObj->setDelGapExt( $del_extn_penalty );


if ( exists $options{'raw'} ) {
  $sEngineObj->setScoreMode( SearchEngineI::basicScoreMode );
}else {
  $sEngineObj->setScoreMode( SearchEngineI::complexityAdjustedScoreMode );
}

# Resolve matrix
my $resolvedMatrix = "";
if ( $options{'matrix'} ) {

  my $matFileName = $options{'matrix'};
  if ( $options{'matrix'} =~ /^(\d+)$/ ) {
    $matFileName = "$1p41g.matrix";
  }elsif ( $options{'matrix'} =~ /^(\d+)p(\d+)$/ ) {
    $matFileName = "$1p$2g.matrix";
  }

  # Simple file reference
  if ( -s $matFileName ) {
    $resolvedMatrix = $matFileName;
  }else {
    # Look through path in priority order
    my @path = ();
    push @path, $ENV{'MATRIX_DIR'} if ( exists $ENV{'MATRIX_DIR'} && -d $options{'MATRIX_DIR'} );
    if ( $engine eq "rmblast" ) {
      push @path, "$FindBin::RealBin/../Matrices/ncbi/nt";
      push @path, $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'} . "/Matrices/ncbi/nt";
    }
    if ( $engine eq "crossmatch" ) {
      push @path, $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'} . "/Matrices/crossmatch";
    }
    foreach my $dir ( @path ) {
      if ( -s "$dir/$matFileName" ) {
        $resolvedMatrix = "$dir/$matFileName";
        last;
      }
    }
  }
  die "ERROR: Could not find matrix ($options{'matrix'})!\n" if ( $resolvedMatrix eq "" );
}else{
  if ( $engine eq "rmblast" ) {
    $resolvedMatrix = "$FindBin::RealBin/../Matrices/ncbi/nt/comparison.matrix";
  }elsif ( $engine eq "crossmatch" ) {
    $resolvedMatrix = "$FindBin::RealBin/../Matrices/crossmatch/comparison.matrix";
  }
}
if ( $resolvedMatrix ne "" ) {
  $sEngineObj->setMatrix( $resolvedMatrix );
}

my $db_seqs;
my $db_size;
if ( $engine eq "rmblast" ) {
  if ( ! -s "$databaseFile.nhr" || $options{'force'} ) {
    if ( $options{'rmb_db_softmask'} ) {
      # The blast database has several hardcoded IDs for a handful of standard filters:
      #    enum EBlast_filter_program {
      #        eBlast_filter_program_not_set      =   0,
      #        eBlast_filter_program_dust         =  10,
      #        eBlast_filter_program_seg          =  20,
      #        eBlast_filter_program_windowmasker =  30,
      #        eBlast_filter_program_repeat       =  40,
      #        eBlast_filter_program_other        = 100,   <<< This is the one we will use
      #        eBlast_filter_program_max          = 255
      #    };
      system(   "$engine_dir/convert2blastmask -in $databaseFile -parse_seqids " 
              . "-masking_algorithm Unknown -masking_options \"all lowercase regions\" "
              . "-outfmt maskinfo_asn1_bin -out $databaseFile.asnb" );
      # Do we need to use blastdb version 4 anymore?
      #      . "-blastdb_version 4 "
      system(   "$engine_dir/makeblastdb -out $databaseFile "
              . "-mask_data $databaseFile.asnb "
              . "-parse_seqids -dbtype nucl -in $databaseFile > "
              . "makedb.log 2>&1" );
    }else {
      # Do we need to use blastdb version 4 anymore?
      #      . "-blastdb_version 4 "
      my $cmd = "$engine_dir/makeblastdb -out $databaseFile "
              . "-parse_seqids -dbtype nucl -in $databaseFile > "
              . "makedb.log 2>&1";
      system( $cmd );
      if ( $? ) {
        printf "\n\nERROR building nucleotide database! makeblastdb command ($cmd) exited with value %d\n\n", $? >> 8;
        system("cat makedb.log");
        exit(1);
      }
    }
  }else {
    unless ( ! $options{'quiet'} ) {
      print "# WARNING: RMBlast database exists for $databaseFile.  Use -force to force rebuilding of the database\n";
    }
  }
  open IN,"$engine_dir/blastdbcmd -db $databaseFile -dbtype nucl -info|" or 
     die "Could not obtain database info using $engine_dir/blastdbcmd";
  while ( <IN> ) {
    # Database: sva_a.fa
    # 	1 sequences; 1,387 total bases
    # 
    # Date: Jun 5, 2023  12:33 PM	Longest sequence: 1,387 bases
    # 
    # BLASTDB Version: 5
    # 
    # Volumes:
    # 	/u1/home/rhubley/projects/cons_thresholds/foo
    if ( /^\s+([\d\,]+)\s+sequences;\s+([\d\,]+)\s+total bases\s*$/ ) {
      $db_seqs = $1;
      $db_size = $2;
    }
  }
  close IN;
  $db_seqs =~ s/,//g;
  $db_size =~ s/,//g;

}elsif ( $engine eq "nhmmer" ) {
  my $dbFile = $databaseFile;
  if ( $options{'fmindex'} ) {
     if ( ! -s "$databaseFile.fm" ) {
       system("$engine_dir/makehmmerdb $databaseFile $databaseFile.fm > makedb.log 2>&1");
     }
     $dbFile = "$databaseFile.fm";
  }
  # If query is a fasta file, treat it as a singlemx search
  if ( $queryFile =~ /.*\.(fa|fna|fasta)/i ) { 
    my $params = " --singlemx ";
    $params .= " --mxfile $resolvedMatrix " if ( $resolvedMatrix );
    # TODO: Add gap probabilities
    $params .= " --cpu $options{'threads'}" if ( $options{'threads'} );
    # TESTING:
    #$params .= " --popen 0.03125 --pextend 0.75 "; # defaults
    # Seems to do better
    #$params .= " --popen 0.021 --pextend 0.52 "; # defaults
    $params .= " $queryFile $dbFile";
    $sEngineObj->setOverrideParameters($params);
  }
}

unless ( $options{'quiet'} ) {
  print "#\n";
  print "# $engine [$engine_ver]\n";
  print "#\n";
  my $cmdLine = $sEngineObj->getParameters();
  $cmdLine =~ s/(.{50}[^\/\s]+)/$1\n\#                  /g;
  if ( $engine eq "rmblast" ) {
    print "#  command_line  : export BLASTMAT=" . dirname($resolvedMatrix) . ";\n";
    print "#                $cmdLine\n";
    print "#  copy_paste    : export BLASTMAT=" . dirname($resolvedMatrix) . "; " . 
          $sEngineObj->getParameters() . "\n"; #if ( $options{'verbose'} );
    print "#  database: $db_seqs sequences, $db_size bp\n";
    print "#  matrix: $resolvedMatrix";
    if ( exists $options{'k_param'} && exists $options{'lambda'} && exists $options{'q_size'} ) { 
      print ", K=$options{'k_param'}, lambda=$options{'lambda'}, q_size=$options{'q_size'}";
    }
    print "\n";
  }else {
    print "#  command_line  : $cmdLine\n";
    print "#  copy_paste    : " . $sEngineObj->getParameters() . "\n";
    #  if ( $options{'verbose'} );
  }
  print "#\n";
}

my $start = time;
my ( $status, $resultCollection ) = $sEngineObj->search();
my $end = time;
unless( $options{'quiet'} ) {
  print "#  search_time(s): " . ( $end - $start + 1 ) . "\n";
  print "#  alignments    : " . $resultCollection->size() . "\n";
  print "#\n";
}

# TODO: add an option to cleanup the databases unless iterative searches are planned

if ( $status )
{
  print STDERR "\nERROR from search engine (", $? >> 8, ") \n";
} else
{
  my $matrixObj;
  if ( $resolvedMatrix ne "" ) {
    $matrixObj = Matrix->new( fileName => $resolvedMatrix );
    # For rescoring purposes we must transpose rmblast-style matrices
    if ( $engine eq "rmblast" ) {
      $matrixObj->transposeMatrix();
    }
  }
  for ( my $k = 0 ; $k < $resultCollection->size() ; $k++ ) {
    my $resultRef = $resultCollection->get( $k );

    #my ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
    #    $resultRef->calcKimuraDivergence( divCpGMod => 1 );
    #
    if ( $engine eq "nhmmer" && $matrixObj && 1 ) {
      # NOTE: This must be a matrix in crossmatch format if matrix is
      #       assymetrical
      my ( $score, $divergence, $cpgsites, $percIns, $percDel,
           $positionScores, $xdrop_fragments, $well_characterized_bases,
           $transisitions, $transversions ) 
                = $resultRef->rescoreAlignment( scoreMatrix => $matrixObj,
                                gapOpenPenalty => $gap_open_penalty,
                                insGapExtensionPenalty => $ins_extn_penalty,
                                delGapExtensionPenalty => $del_extn_penalty,
                                complexityAdjust => 1 );
      if ( $options{'minscore'} ) {
        next if ( $score < $options{'minscore'} );
      }
      $resultRef->setScore($score);
    
    }

    my $bitScore = 0;
    my $e_value = -1;
    my $rawScore;
    if ( ($engine eq "rmblast" || $engine eq "crossmatch") && $options{'k_param'} && $options{'lambda'} && $options{'q_size'} ) {
      $rawScore = $resultRef->getScore();
      my $q_size = $options{'q_size'};
      if ( ! $options{'raw'} ) {
        ## NOTE: For this to work the matrix must be in the orientation expected
        #        by crossmatch.  E.g. not the orientation used by rmblast ( see above where matrixObj
        #        is loaded.
        my ( $raw_score, $divergence, $cpgsites, $percIns, $percDel,
             $positionScores, $xdrop_fragments, $well_characterized_bases,
             $transisitions, $transversions )
                  = $resultRef->rescoreAlignment( scoreMatrix => $matrixObj,
                                  gapOpenPenalty => $gap_open_penalty,
                                  insGapExtensionPenalty => $ins_extn_penalty,
                                  delGapExtensionPenalty => $del_extn_penalty 
                                  );
        $rawScore = $raw_score;
      }
      # TODO: $db_size isn't currently pulled from crossmatch output
      # TODO: We also don't know query size from rmblast
      #
      # NOTE: The search space used in this evalue computation assumes that 
      #       both DNA strands are searched ( ie. search space = 2 * m * n, where
      #       m/n are the single-strand lengths of the subject and query sequences
      #       respectively ).  No attempt is made here to adjust for edge effects
      #       in this calculation. TODO: Consider edge correcting this, however
      #       we would also need alpha/beta in order to do so.
      #
      #       Noteworthy observations: As of NCBI Blast 2.14.1 the evalue calculation
      #       used in blastn doesn't modify the search space when the strand option is
      #       used.  It also appears as if the search space is always assumed to be m*n
      #       , the single stranded search space.  
      #
      #       LAST does appear to adjust the search space for stranded searches.
      #
      $bitScore = sprintf("%0.2f",($options{'lambda'} * $rawScore - log($options{'k_param'})) / log(2));
      $e_value = 2 * $db_size * $q_size * $options{'k_param'} * exp(-$options{'lambda'} * $rawScore);
    }
    
    if ( $options{'alignments'} ) {
      print "" . $resultRef->toStringFormatted( SearchResult::AlignWithQuerySeq );
      if ( $e_value >= 0 ) {
        print "e_value = $e_value\n\n";
      }
    }elsif ( $options{'caf'} ) { 
      # TODO: Move this modified CAF to SearchResult.pm
      #my $matrix = basename($resolvedMatrix);
      my $matrix = $resultRef->getMatrixName();
      print "" . $resultRef->toStringFormatted( SearchResult::CompressedAlignCSV ) . "," . $matrix;
      # Special case for threshold study....append raw score and evalue to caf format
      if ( $options{'lambda'} ) {
        print "," . sprintf("%0.4e", $e_value) . ",$rawScore,$bitScore";
      }
      print "\n";
    }elsif ( $options{'bed3'} ) {
      print "" . $resultRef->getQueryName() . "\t" 
               . ($resultRef->getQueryStart() - 1) . "\t"
               . $resultRef->getQueryEnd(). "\n";
    }elsif ( $options{'bed6'} ) {
      my $orient = "+";
      $orient = "-" if ( $resultRef->getOrientation() );
      print "" . $resultRef->getQueryName() . "\t" 
               . ($resultRef->getQueryStart() - 1) . "\t"
               . $resultRef->getQueryEnd() . "\t" 
               . $resultRef->getSubjName() . "\t" 
               . $resultRef->getScore() . "\t" 
               . $orient . "\n";
    }else {
      print "" . $resultRef->toStringFormatted( SearchResult::OutFileFormat );
      if ( $e_value >= 0 ) {
        print "e_value = $e_value\n";
      }
    }
  }
}


##### NOT PORTED BELOW HERE ####
#die "$USAGE\n$ARGV[0] is not a readable file\n" unless -s $ARGV[0] || $opt_zip;
#my $fastafile1 = shift if $ARGV[0];
#my $fastafile2 = "";
#my $selfcomparison;
#my $subfilecnt;
#my $tempsubjfile;
#while ($ARGV[0]) {
#  if (-s $ARGV[0]) {
#    if ($subfilecnt) {
#      my $extrafile = shift;
#      my $efname = $extrafile;
#      $efname =~ s/^\S+\///; # remove any directory structure
#      if ($tempsubjfile) {
#	$tempsubjfile .= "_$efname";         
#      } else {
#	$tempsubjfile = $fastafile2;
#	$tempsubjfile =~ s/^\S+\///; # remove any directory structure
#	$tempsubjfile = "tempcombinedfiles_$tempsubjfile\_$efname";
#      }
#      system "cat $fastafile2 $extrafile > $tempsubjfile";
#      $fastafile2 = $tempsubjfile;
#    } else {
#      $fastafile2 = shift;
#    }
#    ++$subfilecnt;
#  } else {
#    die "$USAGE\n$ARGV[0] is not a readable file\n"
#  }
#}
#if ($opt_blast && !$fastafile2) {
#  # when cross_match notices that there is not subject file, it assumes
#  # masklevel 101 (it literally sets masklevel to 101) and does not
#  # show self matches (which it does show in "cross_match file file")
#  # Logical as there is no use for an output of self matches only.
#  # rmblast does not work without a subject file. Since this is a
#  # common use, cross_match's behavior is copied by:
#  $fastafile2 = $fastafile1;
#  $opt_level = 101;
#  $selfcomparison = 1;
#  # currently the outputs are different in that with rmblast you get the self matches
#  # while cross_match does not report these
#}
#
#my $matrix;
## default matrices
#if ($opt_blast) {
#  $matrix = "/u1/home/asmit/Matrices/rmblast/nt/ctools.matrix";
##  $matrix = "/u1/local/RepeatModeler-git/Matrices/ncbi/nt/ctools.matrix"; 
## doesn't work for me yet:
##  $matrix = "$FindBin::RealBin/../Matrices/ncbi/nt/ctools.matrix";
#} else {
## for development:
#  $matrix = "/u1/home/asmit/Matrices/ctools.matrix"; # same as /u1/home/asmit/Matrices/matrix
## not precisely the same as /u1/local/RepeatModeler-git/Matrices/crossmatch/ctools.matrix
## scores of R vs A/G/R and Y vs C/T/Y are not symmertrical in ctools.matrix
##  $matrix = "$FindBin::RealBin/../Matrices/crossmatch/ctools.matrix"; 
#}
## perhaps $FindBin::RealBin can find the directory instead of the file. That would be more convenient for development.
#my $matdir = $matrix;
#$matdir =~ s/ctools.matrix$//;
#
#my $command; # will contain the cross_match or rmblastn commandline
#if ($opt_zip) {
#  die "$USAGE\n$opt_zip is not a readable file\n" unless -s $opt_zip;
#  my $sawxmline;
#  open (TEST, $opt_zip);
#  while (<TEST>) {
#    if (/^cross_match/ && !$sawxmline) {
#      unless (/alignments/) {
#	die "With -cg or -perbasescore and -zip you need a cross_match output with alignments as input.\n" if $opt_cg || $opt_perbasescore;
#      }
#      $command = $_;
#      if ($opt_score) {
#	$command =~ s/-score \d+/-score $opt_score/ || $command =~ s/\n$/-score $opt_score\n/;
#      }
#      # extract the matrix used; other parameters do not affect recalculations when changing -p, -cg or -f options
#      if (/-matrix\s+(\S+)/) {
#	$matrix = $1;
#      } else {
#	$matrix = "/u1/home/asmit/Matrices/defaultcrossmatchmatrix";
#      }
#      ++$sawxmline;
#    } elsif ( /^Query file\(s\):\s+(\S+)/) {
#      # not sure yet if the query and subject files are needed for anything 
#      $fastafile1 = $1;
#    } elsif (/^Subject file\(s\):\s+(\S+)/) { # not always there
#      $fastafile2 = $1;
#      last;
#    } elsif (/^Presumed sequence/) {
#      last;
#    } elsif ( $. > 10) {
#      last;
#    }
#  }
#  die "With -zip use an unmodified cross_match output or the .rawxmout output of Xmatch.pl\n" unless $sawxmline;
#  close TEST;
#}
#
#if ( $opt_perbasescore ) {
#  $opt_cg = 1;
#} elsif ($opt_cg) {
#  die "The option -cg has only implemented within the option -perbasescore for now.\n";
#}
#
#my $alignments = "";
#$alignments = "-alignments" if $opt_alignments || $opt_cg;
#
#my $bandwidth = ""; 
#if (defined $opt_bandwidth) {
#  if ($opt_bandwidth < 0) {
#    die "$USAGE\nbandwidth (-b) needs to be positive or 0 for ungapped alignments\n"; 
#  } else {
#    $bandwidth = "-bandwidth $opt_bandwidth";
#  }
#}
#my $gapext = "-gap_ext -5";
#if (defined $opt_del_gap_ext) {
#  die "$USAGE\nEither choose gap_ext (-e) or del_gap_ext (-d), not both" if $opt_extension;
#  $opt_del_gap_ext = -$opt_del_gap_ext if $opt_del_gap_ext > 0;
#  if ($opt_del_gap_ext == 0 ) {
#    die "$USAGE\ndel_gap_ext (-d) should be negative\n";
#  } else {
#    my $ins = $opt_del_gap_ext - 1;
#    $gapext = "-del_gap_ext $opt_del_gap_ext -ins_gap_ext $ins";
#  }
#}
#if (defined $opt_extension) {
#  $opt_extension = -$opt_extension if $opt_extension > 0;
#  if ($opt_extension == 0) {
#    die "$USAGE\ngap_ext (-e) should be negative";
#  } else {
#    $gapext = "-gap_ext $opt_extension";
#  }
#}
#my $gapinit = "-gap_init -25";
#if (defined $opt_gap_init) {
#  $opt_gap_init = -$opt_gap_init if $opt_gap_init > 0;
#  if ($opt_gap_init == 0 ) {
#    die "$USAGE\ngap_init (-g) should be negative \n";
#  } else {
#    $gapinit = "-gap_init $opt_gap_init";
#  }
#}
#
#my $masklevel = ""; # cross_match defaults to -masklevel 80
#if (defined $opt_level) {
#  if ($opt_level < 0  && $opt_level > 101) {
#    die "$USAGE\nmasklevel (-l) is outside the range of 0-101";
#  } else {
#    $masklevel = "-masklevel $opt_level";
#  }
#}
#
#unless ($opt_zip && $matrix) {
#  if ($opt_matrix) {
#    if ($opt_matrix =~ /^\d+$/) {
#      if ($opt_matrix == 14 || $opt_matrix == 18 || $opt_matrix == 20 || $opt_matrix == 25) {
#	$matrix = "$matdir$opt_matrix"."p41g.matrix";
#      } else {
#	die "Can only choose from 14, 18, 20 or 25 for matrix divergence levels\n";
#      }
#    } elsif ($opt_matrix =~ /^\d{2}p\d{2}g?$/) {
#      $matrix = "$matdir$opt_matrix"."g.matrix";
#      $matrix =~ s/gg\.matrix/g\.matrix/;
#      if (-s $matrix) {
#	$matrix = "$matrix";
#      } else {
#	die "$matrix does not seem to exist at $matdir\n";
#      }
#    } elsif (-s $opt_matrix) {
#      $matrix = "$opt_matrix";
#      if ($opt_blast) {
#	$matrix =~ s/\S+\//$matdir/;
#      }
#    } else {
#      $matrix = "$matdir$opt_matrix";
#      if (-s $matrix) {
#	$matrix = "$matrix";
#      } else {
#	die "Cannot find $opt_matrix in current directory or at $matdir\n";
#      }
#    }
#  } elsif ($opt_perbasescore) {
###### CHANGE ONCE PERBASEMATRIX IS IN RMBLAST's DIRECTORY #####
##    $matdir = "/u1/home/asmit/Matrices/rmblast/nt/" if $opt_blast;
#    $matrix = "$matdir"."symmetricNscores0matrix";
#  }
#}
#
## read in matrix
#my @base;
#my %matscore;
#my $j = 0;
#open (INM, "$matrix") || die "File $matrix does not exist\n";
#while (<INM>) {
#  next if /^\#?\s*FREQS/;
#  if (/^\s*[A-Za-z]\s+[A-Za-z]/) {
#    tr/a-z/A-Z/;
#    @base = split;
#  } elsif (/^([A-Za-z])?\s*\-?\d/) {
#    my @score = split;
#    shift @score if $score[0] =~ /[A-Za-z]/;
#    for (my $i = 0; $i <= $#base; ++$i) {
#      # tested with 14p41g.matrix : A in top sequence $i (genomic
#      # copy) aligned to G in bottom sequence $j (e.g. consensus) gets
#      # a penalty of -7 , vice versa -11
#      $matscore{$base[$i]}{$base[$j]} = $score[$i];
#    }
#    ++$j;
#  }
#}
#
#if ($opt_perbasescore) {
#  # establish score of base to self if not present in matrix (cross_match treats them like N)
#  foreach my $b ('A'..'Z') {
#    $matscore{$b}{$b} = $matscore{'N'}{'N'} unless defined $matscore{$b}{$b};
#  }
#  # test for symmetry 
#  my $wtext = "\nThis may lead to unexpected output with the -perbasescore option.
#You could symmetrise your matrix or use the default matrix.\n";
#  for (my $i = 0; $i <= $#base; ++$i) {
#    for (my $j = 0; $j <= $#base; ++$j) {
#      if ( $matscore{$base[$i]}{$base[$j]} != $matscore{$base[$j]}{$base[$i]} ) {
#	die "$matrix is not completely symmetric (seq1 vs seq2 could score
#differently then seq2 vs seq1). For example:
#$base[$i] to $base[$j] scores $matscore{$base[$i]}{$base[$j]} and $base[$j] vs $base[$i] scores $matscore{$base[$j]}{$base[$i]}$wtext";
#      }
#    }
#  } 
#  if ( $matscore{'G'}{'G'} != $matscore{'C'}{'C'} ) {
#    die "$matrix appears strand specific. 
#It scores G:G different ($matscore{'G'}{'G'}) than C:C ($matscore{'C'}{'C'})$wtext";
#  } elsif ( $matscore{'A'}{'A'} != $matscore{'T'}{'T'} ) {
#    die "$matrix appears strand specific.
#It scores A:A different ($matscore{'A'}{'A'}) than T:T ($matscore{'T'}{'T'})$wtext";
#  } 
#}
#$matrix = "-matrix $matrix"; 
#my $CtoTpenalty = $matscore{'C'}{'C'} - $matscore{'C'}{'T'}; # used in calccgadjust sub over and over
#
#my $minmatch = "";
#if ($opt_minmatch) {
#  if ($opt_minmatch < 1) {
#    die "$USAGE\nminmatch (-m) needs to be a positive number";
#  } else {
#    $minmatch = "-minmatch $opt_minmatch";
#  }
#}
#
#my $minscore = "-minscore 200";
#if ($opt_score) {
#  if ($opt_score < 1) {
#    die "$USAGE\nminscore (-s) needs to be a positive number\n";
#  } else {
#    $minscore = "-minscore $opt_score";
#  }
#}
#
## skippable with $opt_zip, but need to rearrange order;
#my $quiet = "";
#$quiet = 1 if $opt_quiet;
#my $raw = "";
#my $wordraw= "";
#my $screen = "";
#if ($opt_blast) {
#  print "\nWARNING: The options -raw, -wordraw or -screen can not yet be used with the rmblastn engine (-blast option).
#The -raw, -wordraw and -screen options are ignored in this run.\n\n" if $opt_raw || $opt_word_raw || $opt_x;
#  print "\nWARNING: Lacking a -raw option, the relative score produced by the -perbasescore option with the rmblastn engine is dependent on the complexity of the sequence. The relative score of low complexity sequences will be lower then they should be.\n\n" if $opt_perbasescore;
#} else {
#  $raw = "-raw" if $opt_raw || $opt_perbasescore;
#  $wordraw = "-word_raw" if $opt_word_raw;
#  $screen = "-screen" if $opt_x;
#}
#
## also skippable with $opt_zip, except for the %ambigadj data; perhaps write that to the rawxmout file?
#
## check for entries with identical names
## count ambiguous bases in query to calculate proper maximum score for $opt_perbasescore
#my %ambigadj;
#my $filecounter;
#my $foundadupname;
#foreach my $bestand ($fastafile1, $fastafile2) {
#  # "bestand" is a Dutch word for file
#  last unless $bestand; # i.e. when $fastafile2 does not exist
#  if ($selfcomparison) {
#    $fastafile2 .= '.modified' if $foundadupname; # $fastafile1 name has been changed to "$fastafile1.modified"
#    last;
#  }
#  my %seenthat;
#  my ($seqname,$seqstr,$allwhatsinbestand);
#  my $modfile = "$bestand.modified";
#  open (IN, "$bestand");
#  while (<IN>) {
#    if ( /^>/ ) {
#      if ($opt_perbasescore &&
#	  !$filecounter && # adjustments are only made for the first file
#	  $seqstr ) {
#	&calcambigadj($seqname,$seqstr);
#	$seqstr = "";
#      }
#      $_ =~ /^>(\S+)/;
#      $seqname = $1;
#      if ($seenthat{$seqname}) {
#	my $newseqname = "$seqname.\[$seenthat{$seqname}\]";
#	my $occ = $seenthat{$seqname} + 1;
#	warn "$seqname occurs $occ times in $bestand. Nr$occ renamed to $newseqname in output and $modfile\n";
#	s/^>$seqname/>$newseqname/;
#	$seqname = $newseqname; # to be fed to calcambigadj
#	$foundadupname = 1;
#      }
#      ++$seenthat{$seqname};
#    } elsif (/\w/) {
#      $seqstr .= $_;
#    }
#    $allwhatsinbestand .= $_;
#  }
#  &calcambigadj($seqname,$seqstr) if $opt_perbasescore && $seqstr;
#  if ($foundadupname) {
#    open (OUTMOD, ">$modfile");
#    print OUTMOD "$allwhatsinbestand";
#    close OUTMOD;
#    if ($bestand eq $fastafile1) {
#      $fastafile1 = $modfile;
#    } elsif ($bestand eq $fastafile2) {
#      $fastafile2 = $modfile;
#    }
#  } 
#  $allwhatsinbestand = "";
#  ++$filecounter;
#}
#
#unless ($opt_zip) {
#  # check if the subject file blast database is up to date
#  if ($opt_blast) {
#    my $ndbfile = "$fastafile2.ndb";
#    if ( -s $ndbfile && (stat($ndbfile))[9] < (stat($fastafile2))[9] ) {
#      unlink ($ndbfile, "$fastafile2.nhr", "$fastafile2.nin", "$fastafile2.nog", "$fastafile2.nos", "$fastafile2.not", "$fastafile2.nsq", "$fastafile2.ntf", "$fastafile2.nto");
#    }
#  }
#
#  # execute cross_match command
#  $command = "$fastafile1 $fastafile2 $alignments $bandwidth $gapext $gapinit $masklevel $matrix $minmatch $minscore $raw $wordraw $screen";
#  if ($opt_blast) {
#    $command = "rmblast.pl $command";
#  } else {
#    $command = "cross_match $command";
#  }
#}
#
#print "$command\n" if $opt_nab;# otherwise already printed to STDOUT by both rmblast.pl and cross_match
#my $rawoutfile = "$fastafile1"."_vsit.rawxmout";
#my $short2 = $fastafile2;
#$short2 =~ s/^\S+\///;
#$rawoutfile = "$fastafile1"."_vs$short2.rawxmout" if $fastafile2;
#if ($opt_zip) {
#  $rawoutfile = $opt_zip;
#} elsif ($quiet) {
#  system "($command > $rawoutfile) >& tempxmatch.stderr";
#} else {
#  system "$command > $rawoutfile";
#}
#
## print sequences in .screen file on one line each
#if ($screen) {
#  my $xdfile = "$fastafile1.screen";
#  if ($opt_zip) {
#    warn "Currently, no (new) masked file $xdfile is created when skipping the cross_match / rmblastn step with the -zip option\n\n";
#  } else {
#    # cross_match prints out a.screen file in 50 letter blocks
#    # something to do with limited memory of computers in the past?
#    # I much prefer sequences to be on one line.
#    my $random = rand(999999);
#    my $tempxdfile = "extraordinarilytemporaryfile$random";
#    open (OUT, ">$tempxdfile");
#    open (IN, "$xdfile");
#    my $linecnt = 0;
#    while (<IN>) {
#      if (/^>/) {
#	if ($linecnt) {
#	  print OUT "\n$_";
#	} else {
#	  print OUT;
#	}
#      } else {
#	chomp;
#	print OUT;
#      }
#      ++$linecnt;
#    }
#    print OUT "\n";
#    close OUT;
#    rename $tempxdfile,$xdfile;
#  }
#}
#
## clean up the temporary file that is a concatenation of subject files
#if ($tempsubjfile) {
#  system "rm $tempsubjfile*";
#}
#
#
## I suppose this is faster than twice reading a file in line by line
#open (IN, "$rawoutfile");
#my @lines = <IN>;
#close IN;
#
#my $print = 1;
#
#if ($opt_perbasescore) {
#  print "The first column is the (raw) score divided by the (raw) score of the query to
#itself times 100. If a non-symmetrical matrix was optioned for the alignment,
#this fraction will depend on the query sequence composition.\n\n"
#}
#
## run through lines to make adjustments with $opt_cg
#my (%addtoscore,%nrmismatches) = ();
#if ($opt_cg) {
#  warn "Making adjustments for CpG sites\n" unless $opt_quiet;
#  my ($xmline,$seq1,$mutline,$seq2,$turn,$leftspace) = ();
#  foreach ( @lines ) {
#    if ( /^\s*\d\d.+\(\d+\)/ ) {
#      if ($xmline) {
#	($addtoscore{$xmline},$nrmismatches{$xmline}) = &calccgadjust($seq1,$mutline,$seq2);
#	$seq1 = "";
#	$mutline = "";
#	$seq2 = "";
#	$turn = 0; #should already be anyway
#      }
#      $xmline = $_;
#    } elsif (/^C?\s+\S+\s+\d+\s([A-Z\-]+)\s\d+/) {
##e.g.
##   AciJub-1.56#LIN         1 ATGGAATATTACTCGGC-ATCAAAAAGAATGAAATCTTGCCATTTGCAAC 49
## or
## C CrasThon-1.62#D       193 TAGGCTGACCAGACGTACCATTTTAAATGGGACAGAAAACAACGGTACTG 144
#      if ($turn) {
#	$seq2 .= $1;
#	$turn = 0;
#      } else {
#	$seq1 .= $1;
#	$turn = 1;
#      }
#      unless ($leftspace) {
#	/^(C?\s*\S+\s+\d+\s+)/;
#	$leftspace = length $1;
#      }
#    } elsif ($turn) {
#      $mutline .= substr($_,$leftspace);
#      chomp $mutline;
#    }
#  }
#  if ( $xmline ) {
#    ($addtoscore{$xmline},$nrmismatches{$xmline}) = &calccgadjust($seq1,$mutline,$seq2);
#  }
#}
#
#my (%max0,%max5,%max6,%max7,%max9,%max10,%max11,%max12);
## to hold the max length of score, leftinquery, subject_name, and leftinsubj
#my @newlines;
#my ($skipselfalignment,$skipforlowscore);
# warn "Adjusting format\n" unless $opt_quiet;
#foreach ( @lines ) { 
#  # the following eliminates the binary-like stuff that appears
#  # when cross_match encounters a base in the matrix it doesn't 
#  # know about (for me usually "Z")
#  if ( /^Num. pairs/ || /^Discrepancy summary/) {
#    $print = 0;
#  } elsif ( !$print && /\(\d+\)/ ) {
#    $print = 1;
#  }
#  if ($print) {
#    if ( /\(\d+\)/ ) {
#      my $cmline = $_;
#      my @bit = split;
#      next unless $bit[11];
#      if ($selfcomparison && $bit[4] eq $bit[8] 
#	  && $bit[5] eq $bit[9] && $bit[7] eq $bit[11]) {
#	# rmblast.pl run of file against itself (no subject file indicated)
#	# cross_match does not report selfmatches, so let's have rmblast.pl also not do this
#	$skipselfalignment = 1;
#      } else {
#	$skipselfalignment = 0;
#	if ( $addtoscore{$_} ) { # only happens wit $opt_cg
#	  # also adjusting the mismatch level in column2/bit[1] by
#	  # considering CpG to TG or CA not as a mismatch each one of
#	  # these represents 10 points in the $addtoscore
#	  $bit[0] += $addtoscore{$_};
#	  s/^(\s*)\d+/$1$bit[0]/;
#	}
#	if ($opt_zip && $opt_score) {
#	  if ($opt_score > $bit[0]) {
#	    $skipforlowscore = 1;
#	  } else {
#	    $skipforlowscore = 0;
#	  }
#	}
#	unless ($opt_original || $bit[1] eq '0.00') {
#	  my $ql = $bit[6] - $bit[5] + 1;
#	  my $mismatches = int ($ql * $bit[1] / 100 + 0.5); 
#	  # proper rounding; if the number was rounded properly to start with
#	  $mismatches = $nrmismatches{$cmline} if defined($nrmismatches{$cmline}); 
#	  # only happens with $opt_cg; can be zero now
#	  my $alignedbases = int ($ql * (1-$bit[3]/100) + 0.5);
#	  my $realmmlevel = sprintf "%4.2f", (100*$mismatches/$alignedbases);
#	  if ($realmmlevel < 10) {
#	    s/^(\s*\d+)\s+$bit[1]/$1  $realmmlevel/;
#	  } else {
#	    s/^(\s*\d+)\s+$bit[1]/$1 $realmmlevel/;
#	  }
#	}
#	# establish maximum width fields for each query name
#	# empty lines separate different queries so it does
#	# not look chaotic to have different widths for each
#        $max0{$bit[4]} = length $bit[0] if !$max0{$bit[4]} || (length $bit[0]) > $max0{$bit[4]};
#        $max5{$bit[4]} = length $bit[5] if !$max5{$bit[4]} || (length $bit[5]) > $max5{$bit[4]};
#        $max6{$bit[4]} = length $bit[6] if !$max6{$bit[4]} || (length $bit[6]) > $max6{$bit[4]};
#        $max7{$bit[4]} = length $bit[7] if !$max7{$bit[4]} || (length $bit[7]) > $max7{$bit[4]};
#	splice @bit,8,0,'+' unless $bit[8] eq 'C' && $bit[12];
#	$max9{$bit[4]} = length $bit[9] if !$max9{$bit[4]} || (length $bit[9]) > $max9{$bit[4]};
#	$max10{$bit[4]} = length $bit[10] if !$max10{$bit[4]} || (length $bit[10]) > $max10{$bit[4]};
#	$max11{$bit[4]} = length $bit[11] if !$max11{$bit[4]} || (length $bit[11]) > $max11{$bit[4]};
#	$max12{$bit[4]} = length $bit[12] if !$max12{$bit[4]} || (length $bit[12]) > $max12{$bit[4]};
#      }
#      push @newlines, $_ unless $skipselfalignment || $skipforlowscore;
#    } else {
#      push @newlines, $_ unless $opt_nab ||$skipselfalignment || $skipforlowscore;
#    }
#  }
#}
#@lines = ();
#my $lastq;
#foreach ( @newlines ) {
#  if ( /\(\d+\)/ ) {
#    # Format space for clearer output. cross_match only adjusts the
#    # distance between query and start position for the longest
#    # length in the column
#    my @bit = split;
#    print "\n" if $opt_nab && (!$lastq || $lastq ne $bit[4]);
#    my ($sp0,$sp5,$sp6,$sp7,$sp9,$sp10,$sp11,$sp12);
#    my $sp1 = "";
#    my $sp2 = "";
#    my $sp3 = "";
#    $sp0 = " " x ($max0{$bit[4]} - (length $bit[0]) );
#    $sp1 = " " if $bit[1] < 10;
#    $sp2 = " " if $bit[2] < 10;
#    $sp3 = " " if $bit[3] < 10;
#    $sp5 = " " x ($max5{$bit[4]}-(length $bit[5]));
#    $sp6 = " " x ($max6{$bit[4]}-(length $bit[6]));
#    $sp7 = " " x ($max7{$bit[4]}-(length $bit[7]));
#    splice @bit,8,0,'+' unless $bit[8] eq 'C' && $bit[12];
#    $sp9 = " " x ($max9{$bit[4]}-(length $bit[9]));
#    $sp10 = " " x ($max10{$bit[4]}-(length $bit[10]));
#    $sp11 = " " x ($max11{$bit[4]}-(length $bit[11]));
#    $sp12 = " " x ($max12{$bit[4]}-(length $bit[12]));
#    $bit[8] = " " if $bit[8] eq '+' && !$opt_nab;
#    $_ = "$sp0$bit[0] $sp1$bit[1] $sp2$bit[2] $sp3$bit[3]  $bit[4]   $sp5$bit[5] $sp6$bit[6] $sp7$bit[7]  $bit[8] $bit[9]$sp9  $sp10$bit[10] $sp11$bit[11] $sp12$bit[12]\n";
#
#    # adding the extra column with per base similarity score
#    if ($opt_perbasescore) {
#      $bit[7] =~ tr/()//d;
#      # non-complexity adjusted score of the full sequence against itself
#      my $maxsc = ($bit[6] + $bit[7])*$matscore{'A'}{'A'};
#      # modify expected score vs itself given ambiguous bases in the sequence
#      $maxsc -= $ambigadj{$bit[4]} if $ambigadj{$bit[4]};
#      my $perbase = sprintf("%4.1f", 100*$bit[0]/$maxsc);
#      if ($perbase >= 100) {
#	print "100.  $_";
#      } else {
#	print "$perbase  $_";
#      }
#    } else {
#      print "$_";
#    }
#    $lastq = $bit[4];
#  } else {
#    print "$_";
#  }
#}
#
#unlink $rawoutfile;
#
#
#
#
#
######################### S U B R O U T I N E S ############################
#
###-------------------------------------------------------------------------##
### Use: my _privateMethod( $parameter => value );
###
###      $parameter       : A parameter to the method
###
###  Returns
###      Methods with the prefix "_" are conventionally considered
###      private.  This bit-o-documentation is not formatted to
###      print out when perldoc is run on this file.
###
###-------------------------------------------------------------------------##
## This appears to come up with a score adjustment
#sub calccgadjust {
#  my $sq1 = shift; # full aligned sequence (in upper case, with dashes and all) of query
#  my $mutationline = shift; # with i v - and ?
#  my $sq2 = shift; # full aligned sequence of target
#  my $adjustment = 0; # has to set to zero
#
#  my $mismatches = ( $mutationline =~ tr/iv// );
#  # Note that this does not count question marks and most aligned
#  # ambiguous bases are not counted as mismatches. This is what I want
#  # in $opt_perbase , but may not be always right when calccgadjust is
#  # used outside -perbase For example, K:A, W:S, M:G etc should always
#  # be counted as a mismatch. (R:Y already gets a "v" by
#  # cross_match). Those are some rare situations, I suppose.
#
#  # are there CpGs in the aligned part of the query?
#  my $sq1cp = $sq1;
#  my $mutlinecp = $mutationline;
#  my ($len1,$len2);
#  while ($sq1cp =~ s/^(\S*?)C(\-*)G//) { # query has a CG pair; change it to cg
#    $len1 = length $1; #distance to C of first remaining CG
#    $len2 = $len1 + (length $2) + 1; # distance to G of that CG (incl the C)
#    # (length $2) representing the number of dashes between C and G
#
#    if ( $mutlinecp =~ /^.{$len1}i/ || $mutlinecp =~ /^.{$len2}i/ ) {  
#      # cross_match indicated a transition; there will be only one adjustment for "ii" (i.e. CG->AT) 
#      $adjustment += 4*$CtoTpenalty/5;
#      --$mismatches;
#      # $CtoTpenalty equals $matscore{'C'}{'C'} - $matscore{'C'}{'T'}
#      # what this does is reduce the original difference of a match
#      # and a mismatch to the total score to 1/5th of that difference
#
#    } elsif ( $mutlinecp =~ /^.{$len1}\?\s/ || $mutlinecp =~ /^.{$len2}\?/ ) {
#      # this is a CpG lined up to NG or CN in the second sequence;
#      # since CpGs have a lot of variability NG or CN are common
#      # outcomes The adjustment is small, because in the -perbase mode
#      # we do not want to make CN and CG calls equally good
#      $adjustment += $CtoTpenalty/3;
#      # no mismatch subtraction as the question marks were not counted
#    # (could do different adjustments when CpG options works outside perbase option)
#    }
#    $mutlinecp =~ s/^.{$len2}.//; #the extra period for the removed G
#  }
#  # are there CpGs in the aligned part of the target? 
#  # doesn't matter if the query has a CpG there too, as no substitutions will be reported
#  my $sq2cp = $sq2;
#  $mutlinecp = $mutationline;
#  while ($sq2 =~ s/^(\S*?)C(\-*)G//) {
#    my $len1 = length $1;
#    my $len2 = $len1 + (length $2) + 1;
#    if ( $mutlinecp =~ /^.{$len1}i/ || $mutlinecp =~ /^.{$len2}i/ )  {
#      --$mismatches;
#      $adjustment += 9*$CtoTpenalty/10;
#      # difference with match is now 1/10 instead of 1/5 above, since we favor 
#      # models with CpG calls over ones with CA and TG calls
#      # no adjustment for NG or CN opposite CG in target, as Ns were
#      # already adjusted (+9 or +10) when counting the ambiguous bases in the query
#    }
#    $mutlinecp =~ s/^.{$len2}.//;
#  }
#  # adjust for CA opposite TG
#  # to speed things up, we skip looking for situations like:
##  yoh            1 GCATTGCATGGTTTAGATAGCTTCGCCATT---GTGTAGTGGTA 41
##                                                i---i          
##  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCTTAATGTAGTGGTA 44
#  # also as it seems less likely that this involved a CpG transition 
#  # It could  if the insert in "yah" starts with an A or ends with a C,
#  # but it appears that cross_match aligns the transitions on one side
#  # of a gap of more than 1 bp:
##  yoh            1 GCATTGCATGGTTTAGATAGCTTCGCCAT---TGTGTAGTGGTA 41
##                                                ---ii          
##  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCTTCATGTAGTGGTA 44
#  # and
##  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCACAATGTAGTGGTA 44
##                                                ii---          
##  yoh            1 GCATTGCATGGTTTAGATAGCTTCGCCATTG---TGTAGTGGTA 41
#  # but:
##  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCAATGTAGTGGTA 42
##                                                i-i          
##  yo             1 GCATTGCATGGTTTAGATAGCTTCGCCATT-GTGTAGTGGTA 41
#
## For CG vs CA/TG 
##  AciJub-6.1674#L         1 GATGGCGGAACAGCATGGAAGTTTTTTGC--GTCTCTCGTCCATGAAATA 48
##                                                       vi--     v             
##  ParHer-4.26#LIN        12 GATGGCGGAACAGCATGGAAGTTTTTTTTCTGTCTCGCGTCCATGAAATA 61
## so I'll keep the check for multi-base gaps above, 
#
#  $sq1cp = $sq1;
#  $mutlinecp = $mutationline;
#  while ($sq1cp =~ s/^(\S*?)C(\-?)A//) { # *? for the shortest distance
#    my $len = length $1;
#    my $gap = $2;
#    if (!$gap && $mutlinecp =~ /^.{$len}ii/ 
#	 || $gap && $mutlinecp =~ /^.{$len}i\-i/ ) {
#      --$mismatches;
#      $adjustment += 9*$CtoTpenalty/10;
#      # same adjustment as for CG opposite CA or TG, since here one transition
#      # represents a CpG change, the other a normal transition
#    }
#    $len += 2; # for the C and A
#    ++$len if $gap;
#    $mutlinecp =~ s/^.{$len}//;
#  }
#  # adjust for TG opposite CA
#
#  $sq1cp = $sq1;
#  $mutlinecp = $mutationline;
#  while ($sq1cp =~ s/^(\S*?)T(\-?)G//) {
#    my $len = length $1;
#    my $gap = $2;
#    if (!$gap && $mutlinecp =~ /^.{$len}ii/
#         || $gap && $mutlinecp =~ /^.{$len}i\-i/ ) {
#      --$mismatches;
#      $adjustment += 9*$CtoTpenalty/10;
#    }
#    $len += 2;
#    ++$len if $gap;
#    $mutlinecp =~ s/^.{$len}//;
#  }
#  $adjustment = int($adjustment + 0.49); # rounding 100.5 down to 100 ; feels better;-)
#  return ($adjustment,$mismatches);
#}
#
#
#
###-------------------------------------------------------------------------##
### Use: my _privateMethod( $parameter => value );
###
###      $parameter       : A parameter to the method
###
###  Returns
###      Methods with the prefix "_" are conventionally considered
###      private.  This bit-o-documentation is not formatted to
###      print out when perldoc is run on this file.
###
###-------------------------------------------------------------------------##
#sub calcambigadj {
#  my $name = shift;
#  my $seq = shift;
#  $seq =~ s/\s//g;
#  $seq =~ tr/a-z/A-Z/;
##  my $gc = ($seq =~ tr/GC//);                                                                        
##  my $at = ($seq =~ tr/AT//);                                                                        
##  my $gclevel = $gc/($at+$gc); # not /length, given ambiguous bases                                  
#  my $match = $matscore{'A'}{'A'};
#  $seq =~ tr/AT//d;
#  my $gc = ($seq =~ tr/GC//d);
#  if ($matscore{'G'}{'G'} != $match) {
#    # earlier Xmatch.pl forbids C:C and G:G or A:A and T:T to have different scores                   
#    $ambigadj{$name} += ($match - $matscore{'G'}{'G'})*$gc;
#  }
#  my $n= ($seq =~ tr/N//d);
#  $ambigadj{$name} += ($match - $matscore{'N'}{'N'})*$n if $n;
#  if ($seq) { # most seqs do not have other ambigs than N                                             
#    my @ambigs = qw (R Y K M S W X B D H V Z);
#    foreach my $a (@ambigs) {
#      $ambigadj{$name} += ($match - $matscore{$a}{$a})*($seq =~ s/$a//g);
#      # cannot use faster tr as it doesnt interpret variables                                         
#      last unless $seq;
#    }
#  }
#  if ($seq) {
#    warn "\nNOTE: $name has irregular bases $seq\n\n";
#  }
#}
#
#

 


