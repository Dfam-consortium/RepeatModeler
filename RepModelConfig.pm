#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) RepModelConfig.pm
##  Author:
##      Arian Smit <asmit@systemsbiology.org>
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      This is the main configuration file for the RepeatModeler
##      program suite.  Before you can run the programs included
##      in this package you will need to run the "./configure" program
##      to modify this file or manually edit in an editor.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2004-2019 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
package RepModelConfig;
use FindBin;
use File::Basename;
use Data::Dumper;
require Exporter;
@EXPORT_OK   = qw();
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );
@ISA         = qw(Exporter);
my $CLASS = "RepModelConfig";

BEGIN {
##----------------------------------------------------------------------##
##     CONFIGURE THE FOLLOWING PARAMETERS FOR YOUR INSTALLATION         ##
##                                                                      ##
##  This file may be hand edited or preferably, modified with the       ##
##  the use of the "configure" script.                                  ##
##                                                                      ##
##  In the following section default values for paths/programs          ##
##  may be hardcoded. Each parameter appears in a block of text         ##
##  below as:                                                           ##
##                                                                      ##
##    PARAMETER_NAME => {                                               ##
##                 ...                                                  ##
##              value => "/usr/local/bin/foo"                           ##
##                      }                                               ##
##                                                                      ##
##  only change the "value" field for each parameter. If you are        ##
##  unsure of how to edit this file, simply use the the ./configure     ##
##  script provided with the package to perform the same task.  For     ##
##  more details on configuring this package please see the README      ##
##  file.                                                               ##
##                                                                      ##
  ## STCFG --do-not-remove--
  $configuration = {
          'DFAM_CURATOR_DIR' => {
                                   'command_line_override' => 'dfam_curator_dir',
                                   'description' => 'The path to the directory holding the dfam-curator programs ( 0.6.6 or higher ).  RepeatModeler builds and extends family consensi with its te-composer, in place of the Refiner script it once ran, and the LTR structural pipeline builds its consensus models and seed alignments with it.',
                                   'environment_override' => 'DFAM_CURATOR_DIR',
                                   'expected_binaries' => [
                                                            'te-composer'
                                                          ],
                                   'expected_files' => [],
                                   'param_type' => 'directory',
                                   'required' => 1,
                                   'value' => '/home/rhubley/projects/Claude/dfam-curator/target/release'
                                 },
          'FAMDB_DIR' => {
                           'command_line_override' => 'famdb_dir',
                           'description' => 'The path to the installation of the FamDB program (V3.0.0) or higher.',
                           'environment_override' => 'FAMDB_DIR',
                           'expected_binaries' => [
                                                    'famdb.py'
                                                  ],
                           'expected_files' => [],
                           'param_type' => 'directory',
                           'required' => 1,
                           'value' => '/usr/local/FamDB-Dfam-4.0_RB'
                         },
          'LTR_FORAGE_DIR' => {
                                'command_line_override' => 'ltr_forage_dir',
                                'description' => 'The path to the directory holding the ltr-forage program ( 0.2.0 or higher ): the bin/ directory of a release, or target/release of a source build.  ltr-forage\'s ltr_library.py is found beside it ( a release\'s bin/ltr_library.py ), in its pipeline/ subdirectory, or in ../../pipeline for target/release.  The LTR structural pipeline ( -LTRStruct ) builds its families with ltr-forage\'s family-library driver.',
                                'environment_override' => 'LTR_FORAGE_DIR',
                                'expected_binaries' => [
                                                         'ltr-forage'
                                                       ],
                                'expected_files' => [],
                                'param_type' => 'directory',
                                'required' => 1,
                                'value' => ''
                              },
          'MMSEQS_DIR' => {
                            'command_line_override' => 'mmseqs_dir',
                            'description' => 'The path to the bin directory of the MMseqs2 installation.  nail uses mmseqs to seed its protein search.  nail will absorb mmseqs in a future release, and this setting goes with it.',
                            'environment_override' => 'MMSEQS_DIR',
                            'expected_binaries' => [
                                                     'mmseqs'
                                                   ],
                            'expected_files' => [],
                            'param_type' => 'directory',
                            'required' => 1,
                            'value' => '/usr/local/MMseqs2-18/bin'
                          },
          'NAIL_DIR' => {
                          'command_line_override' => 'nail_dir',
                          'description' => 'The path to the directory holding the nail program.  RepeatClassifier uses nail for its protein search, in place of the NCBI blastx it once called.',
                          'environment_override' => 'NAIL_DIR',
                          'expected_binaries' => [
                                                   'nail'
                                                 ],
                          'expected_files' => [],
                          'param_type' => 'directory',
                          'required' => 1,
                          'value' => '/usr/local/nail-v0.7.1'
                        },
          'RECON_DIR' => {
                           'command_line_override' => 'recon_dir',
                           'description' => 'The path to the installation of the RECON de-novo repeatfinding program.',
                           'environment_override' => 'RECON_DIR',
                           'expected_binaries' => [
                                                    'eledef',
                                                    'eleredef'
                                                  ],
                           'expected_files' => [],
                           'param_type' => 'directory',
                           'required' => 1,
                           'value' => '/usr/local/RECON'
                         },
          'REPEATAFTERME_DIR' => {
                                   'command_line_override' => 'repeatafterme_dir',
                                   'description' => 'The path to the installation of the RepeatAfterMe ( 0.0.7 or higher ) extension program.',
                                   'environment_override' => 'REPEATAFTERME_DIR',
                                   'expected_binaries' => [
                                                            'RAMExtend'
                                                          ],
                                   'expected_files' => [],
                                   'param_type' => 'directory',
                                   'required' => 1,
                                   'value' => '/usr/local/RepeatAfterMe'
                                 },
          'REPEATMASKER_DIR' => {
                                  'command_line_override' => 'repeatmasker_dir',
                                  'description' => 'The path to the installation of RepeatMasker (4.2.5 or higher).',
                                  'environment_override' => 'REPEATMASKER_DIR',
                                  'expected_binaries' => [
                                                           'RepeatMasker'
                                                         ],
                                  'expected_files' => [],
                                  'param_type' => 'directory',
                                  'required' => 1,
                                  'value' => '/home/rhubley/projects/Claude/RepeatMasker'
                                },
          'RMBLAST_DIR' => {
                             'command_line_override' => 'rmblast_dir',
                             'description' => 'The path to the bin directory of RMBlast (2.17.1 or higher, including the 3.x series).',
                             'environment_override' => 'RMBLAST_DIR',
                             'expected_binaries' => [
                                                      'rmblastn',
                                                      'dustmasker'
                                                    ],
                             'expected_files' => [],
                             'param_type' => 'directory',
                             'required' => 1,
                             'value' => '/usr/local/rmblast-3.0.9/bin'
                           },
          'RSCOUT_DIR' => {
                            'command_line_override' => 'rscout_dir',
                            'description' => 'The path to the installation of the RepeatScout ( 1.0.7 or higher ) de-novo repeatfinding program.',
                            'environment_override' => 'RSCOUT_DIR',
                            'expected_binaries' => [
                                                     'RepeatScout',
                                                     'build_lmer_table'
                                                   ],
                            'expected_files' => [],
                            'param_type' => 'directory',
                            'required' => 1,
                            'value' => '/home/rhubley/projects/RepeatScout-1.0.7'
                          },
          'SUFR_DIR' => {
                          'command_line_override' => 'sufr_dir',
                          'description' => 'The path to the directory holding the sufr suffix-array program ( 0.8.0 or higher ).  The LTR structural pipeline ( -LTRStruct ) indexes the genome with it for ltr-forage.',
                          'environment_override' => 'SUFR_DIR',
                          'expected_binaries' => [
                                                   'sufr'
                                                 ],
                          'expected_files' => [],
                          'param_type' => 'directory',
                          'required' => 1,
                          'value' => ''
                        },
          'TRF_DIR' => {
                         'command_line_override' => 'trf_dir',
                         'description' => 'The full path to TRF program.  TRF must be named "trf". ( 4.0.9 or higher )',
                         'environment_override' => 'TRF_DIR',
                         'expected_binaries' => [
                                                  'trf'
                                                ],
                         'expected_files' => [],
                         'param_type' => 'directory',
                         'required' => 1,
                         'value' => '/usr/local/bin'
                       },
          'UCSCTOOLS_DIR' => {
                               'command_line_override' => 'ucsctools_dir',
                               'description' => 'The path to the installation directory of the UCSC TwoBit Tools (twoBitToFa, faToTwoBit, twoBitInfo etc).',
                               'environment_override' => 'UCSCTOOLS_PRGM',
                               'expected_binaries' => [
                                                        'twoBitToFa',
                                                        'faToTwoBit',
                                                        'twoBitInfo'
                                                      ],
                               'expected_files' => [],
                               'param_type' => 'directory',
                               'required' => 1,
                               'value' => '/usr/local/ucscTools'
                             }
        };

  ## EDCFG --do-not-remove--

##                                                                      ##
##                      END CONFIGURATION AREA                          ##
##----------------------------------------------------------------------##
##----------------------------------------------------------------------##
##  Do not edit below this line                                         ##
##----------------------------------------------------------------------##

  #
  # Apply environment overrides now, while still inside BEGIN, so that
  # "use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}..."
  # in the scripts sees them.  resolveConfiguration() applies them again
  # at run time together with the command line overrides.
  #
  foreach my $param ( keys %$configuration ) {
    my $envName = $configuration->{$param}->{'environment_override'};
    $configuration->{$param}->{'value'} = $ENV{$envName}
        if ( defined $envName && exists $ENV{$envName} );
  }

  #
  # Current version of the software
  #
  $VERSION = "2.0.9";

  #
  # Set this flag to default to debug mode for the entire package
  #
  $DEBUGALL = 0;

  #
  # Prompt for a specific parameter and update the object
  #
  sub promptForParam {
    my $param     = shift;
    my $screenHdr = shift;

    if ( !exists $configuration->{$param} ) {
      return;
    }

    # Grab defaults
    my $defaultValue;
    if ( $configuration->{$param}->{'param_type'} eq "directory" ) {
      if ( exists $configuration->{$param}->{'expected_binaries'}
           && @{ $configuration->{$param}->{'expected_binaries'} } )
      {
        my $binary = $configuration->{$param}->{'expected_binaries'}->[ 0 ];
        $defaultValue = `/usr/bin/which $binary`;
        if ( $defaultValue !~ /\/usr\/bin\/which:/ ) {
          $defaultValue =~ s/[\n\r\s]+//g;
          $defaultValue =~ s/^(.*)\/$binary/$1/;
        }
        else {
          $defaultValue = "";
        }
      }
      if (    $defaultValue eq ""
           && exists $configuration->{$param}->{'value'}
           && -d $configuration->{$param}->{'value'} )
      {
        $defaultValue = $configuration->{$param}->{'value'};
      }
    }
    elsif ( $configuration->{$param}->{'param_type'} eq "program" ) {

      # The program type is used in cases where a single
      # script/binary is referenced and may not have the
      # exact name we expect.  TRF is a good example of this
      # as the binary is often distributed with names like:
      # trf409.linux64 etc..
      if ( exists $configuration->{$param}->{'value'} ) {
        my ( $binary, $dirs, $suffix ) =
            fileparse( $configuration->{$param}->{'value'} );
        $defaultValue = `/usr/bin/which $binary`;
        if ( $defaultValue !~ /\/usr\/bin\/which:/ ) {
          $defaultValue =~ s/[\n\r\s]+//g;
        }
        else {
          $defaultValue = $configuration->{$param}->{'value'};
        }
      }
    }

    my $value = "";
    my $validParam;
    do {
      $validParam = 1;
      system( "clear" );
      if ( $screenHdr ) {
        print "$screenHdr\n";
      }
      else {
        print "\n\n\n\n";
      }
      print "" . $configuration->{$param}->{'description'} . "\n";

      # Prompt and get the value
      if ( $defaultValue ) {
        print "$param [$defaultValue]: ";
      }
      else {
        print "$param: ";
      }
      $value = <STDIN>;
      $value =~ s/[\n\r]+//g;
      if ( $value eq "" && $defaultValue ) {
        $value = $defaultValue;
      }

      if ( $configuration->{$param}->{'param_type'} eq "directory" ) {
        if ( -d $value ) {
          foreach my $file (
                          @{ $configuration->{$param}->{'expected_files'} // [] } )
          {
            if ( !-s "$value/$file" ) {
              print "\nCould not find the required file \"$file\" inside\n"
                  . "the directory \"$value\"!\n\n";
              $validParam = 0;
              last;
            }
          }
          foreach my $binary (
                          @{ $configuration->{$param}->{'expected_binaries'} // [] } )
          {
            if ( !-x "$value/$binary" ) {
              print "\nCould not find the required program \"$binary\" inside\n"
                  . "the directory \"$value\"!\n\n";
              $validParam = 0;
              last;
            }
            elsif ( -d "$value/$binary" ) {
              print "\nCould not find the required program \"$binary\" inside\n"
                  . "the directory \"$value\"!  It appears to be the name of a\n"
                  . "subdirectory.\n\n";
              $validParam = 0;
              last;
            }
          }
        }
        else {
          print "\nCould not find the \"$value\" directory.\n\n";
          $validParam = 0;
        }
      }
      elsif ( $configuration->{$param}->{'param_type'} eq "program" ) {
        if ( !-x $value ) {
          print "\nThe program \"$value\" doesn't appear to exist\n"
              . "or it's not executable!\n\n";
          $validParam = 0;
        }
        elsif ( -d $value ) {
          print "\nThe value \"$value\" appears to be a directory rather\n"
              . "than an executable binary or script!\n\n";
          $validParam = 0;
        }
      }

      if ( $validParam == 0 ) {
        print "<PRESS ENTER TO CONTINUE, CTRL-C TO BREAK>\n";
        <STDIN>;
      }
    } while ( $validParam == 0 );
    $configuration->{$param}->{'value'} = $value;
    #my $version = &getDependencyVersion($param);
    #if ( $version ) {
    #  $configuration->{$param}->{'version'} = $version;
    #}
  }

  #
  # Validate parameter
  #
  sub validateParam {
    my $param       = shift;
    my $new_setting = shift;

    if ( !exists $configuration->{$param} ) {
      return 0;
    }

    my $value = $configuration->{$param}->{'value'};
    $value = $new_setting if ( defined $new_setting );

    # Always assume the "good" in parameters...
    my $validParam = 1;
    if ( $configuration->{$param}->{'param_type'} eq "directory" ) {
      if ( -d $value ) {
        foreach my $file (
                     @{ $configuration->{$param}->{'expected_files'} // [] } )
        {
          if ( !-s "$value/$file" ) {
            $validParam = 0;
            last;
          }
        }
        foreach my $binary (
                     @{ $configuration->{$param}->{'expected_binaries'} // [] } )
        {
          if ( !-x "$value/$binary" ) {
            $validParam = 0;
            last;
          }
          elsif ( -d "$value/$binary" ) {
            $validParam = 0;
            last;
          }
        }
      }
      else {
        $validParam = 0;
      }
    }
    elsif ( $configuration->{$param}->{'param_type'} eq "program" ) {
      if ( !-x $value ) {
        $validParam = 0;
      }
      elsif ( -d $value ) {
        $validParam = 0;
      }
    }

    # TODO Validate versions
    #if ( $validParam ) {
    #  my $version = &getDependencyVersion($param);
    #  if ( $version ) {
    #    $configuration->{$param}->{'version'} = $version;
    #  }
    #}

    return $validParam;
  }

  #
  # Update this file ( beware: self modifying code ) new
  # paramter settings.
  #
  sub updateConfigFile {
    open IN,  "<$CLASS.pm"     or die;
    open OUT, ">new-$CLASS.pm" or die;
    my $inCfg;
    $Data::Dumper::Sortkeys = 1;
    while ( <IN> ) {
      if ( /##\s+STCFG/ ) {
        $inCfg = 1;
        print OUT;
        my $cStr = Dumper( $configuration );
        $cStr =~ s/\$VAR1/  \$configuration/;
        print OUT "$cStr\n";
      }
      elsif ( /##\s+EDCFG/ ) {
        $inCfg = 0;
        print OUT;
      }
      elsif ( !$inCfg ) {
        print OUT;
      }
    }
    close IN;
    close OUT;
    rename( "$CLASS.pm",     "$CLASS.pm.bak" );
    rename( "new-$CLASS.pm", "$CLASS.pm" );
  }

  #
  # Create a GetOpt list for the command-line parameters defined
  # in this configuration file.  These may be appended to a program's
  # existing GetOpt parameter as:
  #
  #     push @getopt_args, RepModelConfig::getCommandLineOptions();
  #
  sub getCommandLineOptions {
    my @options = ();
    foreach my $param ( keys %$configuration ) {
      if ( exists $configuration->{$param}->{'command_line_override'} ) {
        push @options,
            "-" . $configuration->{$param}->{'command_line_override'} . "=s";
      }
    }
    return @options;
  }

  #
  # Get POD documentation to add to the existing program POD stored in
  # the main script.
  #
  sub getPOD {
    my $pod_str;
    foreach my $param ( keys %$configuration ) {
      if ( exists $configuration->{$param}->{'command_line_override'} ) {
        $pod_str .= "=item -"
            . $configuration->{$param}->{'command_line_override'}
            . " <string>\n\n";
        $pod_str .= $configuration->{$param}->{'description'} . "\n\n";
      }
    }
    if ( $pod_str ) {
      return ( "\n=over 4\n\n" . $pod_str . "=back\n\n" );
    }
    return;
  }

  #
  # After GetOpt has filled in the options hash simply pass it to
  # this function to perform resolution.  The following precedence
  # is used:
  #
  #    1. Command Line Parameter
  #    2. Environment Variable
  #    3. Configuration File
  #
  # This will update the configuration{param}->{'value'} for use
  # in the main program.
  #
  sub resolveConfiguration {
    my $options = shift;

    foreach my $param ( keys %$configuration ) {
      if (
        exists $options->{ $configuration->{$param}->{'command_line_override'} }
          )
      {
        $configuration->{$param}->{'value'} =
            $options->{ $configuration->{$param}->{'command_line_override'} };
      }
      elsif (
             exists $ENV{ $configuration->{$param}->{'environment_override'} } )
      {
        $configuration->{$param}->{'value'} =
            $ENV{ $configuration->{$param}->{'environment_override'} };
      }
    }
  }

  #
  # Locate a program from the NCBI BLAST+ toolkit ( only blastdbcmd, to
  # convert a database an earlier release built ).  A 2.x RMBlast ships the
  # toolkit; the 3.x series does not.
  #
  #   my $path = RepModelConfig::findNCBIProgram( "blastdbcmd" );
  #
  # Returns the path of the first executable found in RMBLAST_DIR, then the
  # directories on PATH, or undef.
  #
  sub findNCBIProgram {
    my $name = shift;

    my @dirs = ();
    foreach my $param ( "RMBLAST_DIR" ) {
      my $dir = $configuration->{$param}->{'value'};
      push @dirs, $dir if ( defined $dir && $dir ne "" );
    }
    push @dirs, split( /:/, $ENV{'PATH'} // "" );

    foreach my $dir ( @dirs ) {
      next if ( $dir eq "" );
      my $path = "$dir/$name";
      return $path if ( -x $path && !-d $path );
    }
    return undef;
  }

  #
  # The message to die with when findNCBIProgram() comes up empty.
  #
  sub ncbiProgramMissingMessage {
    my $name    = shift;
    my $purpose = shift;

    return "Could not find the NCBI BLAST+ program '$name', which $purpose.\n"
        . "It was not found in RMBLAST_DIR or on PATH.  Install the NCBI\n"
        . "BLAST+ toolkit ( or a 2.x RMBlast ) and put its bin directory on\n"
        . "PATH.\n";
  }

  #
  # Locate a program RepeatModeler does not configure a directory for, such
  # as cross_match or nhmmer.  Searches PATH and returns the full path, or
  # undef if it finds nothing executable.
  #
  sub findProgramOnPath {
    my $name = shift;

    foreach my $dir ( split( /:/, $ENV{'PATH'} // "" ) ) {
      next if ( $dir eq "" );
      my $path = "$dir/$name";
      return $path if ( -x $path && !-d $path );
    }
    return undef;
  }

  #
  # The message to die with when findProgramOnPath() comes up empty.
  #
  sub programMissingMessage {
    my $name    = shift;
    my $purpose = shift;

    return "Could not find the program '$name', which $purpose.\n"
        . "It is not on your PATH.  Install it and add its directory to\n"
        . "your PATH.\n";
  }

  #
  # Locate the mmseqs binary nail uses for seeding: MMSEQS_DIR first
  # ( config file, MMSEQS_DIR environment variable or -mmseqs_dir ),
  # then PATH.  Returns the full path, or undef.
  #
  sub findMMseqs {
    my $dir = $configuration->{'MMSEQS_DIR'}->{'value'};
    if ( defined $dir && $dir ne "" ) {
      my $path = "$dir/mmseqs";
      return $path if ( -x $path && !-d $path );
    }
    return findProgramOnPath( "mmseqs" );
  }

  #
  # Check that an mmseqs binary runs on this machine.  The MMseqs2 cmake
  # build uses -march=native by default, and mmseqs does not check the CPU
  # at run time, so a copy built on a newer CPU dies with SIGILL on an
  # older one.  main() has no AVX instructions, so "mmseqs version" can
  # succeed on a CPU the build does not support.  This runs a one-sequence
  # search instead, which goes through the prefilter and alignment where
  # the SIMD code is.  It takes about three seconds.
  #
  # Returns "" if the search succeeds, otherwise a message explaining the
  # failure.
  #
  sub checkMMseqs {
    my $mmseqs = shift;

    require File::Temp;
    my $tmpDir = File::Temp::tempdir( "mmseqsCheck-XXXXXX", TMPDIR => 1,
                                      CLEANUP => 1 );
    my $query = "$tmpDir/query.fa";
    open my $FA, ">", $query or return "Could not write $query: $!\n";
    print $FA ">test\n"
        . "MKVLAAGIVGLLLAQSAFAEETPKWRCSHHHNLLQEGVRDLCRKAGLYGDPVEYSIVFV\n";
    close $FA;

    # Fork and exec rather than go through a shell, so a SIGILL shows up
    # in $? as the signal.  The child works inside $tmpDir so any core
    # file is removed with it.
    my $pid = fork();
    return "Could not fork to run $mmseqs: $!\n" unless ( defined $pid );
    if ( $pid == 0 ) {
      chdir $tmpDir;
      open STDOUT, ">", "$tmpDir/check.log";
      open STDERR, ">&", \*STDOUT;
      exec $mmseqs, "easy-search", $query, $query, "$tmpDir/out.m8",
          "$tmpDir/tmp", "--threads", "1";
      exit 127;
    }
    waitpid( $pid, 0 );
    my $status = $?;
    return "" if ( $status == 0 && -e "$tmpDir/out.m8" );

    # If a wrapper script runs mmseqs without exec, the shell reports the
    # SIGILL as exit code 132 ( 128 + 4 ).
    my $signal = $status & 127;
    my $exit   = $status >> 8;
    if ( $signal == 4 || $exit == 132 ) {
      return "$mmseqs died with an illegal instruction ( SIGILL ).  It\n"
          . "was built for CPU features this machine lacks, most likely AVX2 or\n"
          . "AVX-512 from a -march=native build on a newer CPU.  Either\n"
          . "install a build that runs on every machine you use ( the\n"
          . "mmseqs-linux-sse41 release, or build with cmake\n"
          . "-DHAVE_SSE4_1=1 ), or install one build per CPU\n"
          . "type behind a wrapper script that picks the right one from\n"
          . "/proc/cpuinfo.\n";
    }
    if ( $exit == 127 ) {
      return "Could not run $mmseqs.\n";
    }

    my $log = "";
    if ( open my $LOG, "<", "$tmpDir/check.log" ) {
      my @lines = <$LOG>;
      close $LOG;
      $log = join( "", @lines > 10 ? @lines[ -10 .. -1 ] : @lines );
    }
    return "A test search with $mmseqs failed ( "
        . ( $signal ? "signal $signal" : "exit status $exit" ) . " ).\n"
        . ( $log ? "Last lines of its output:\n$log" : "" );
  }

  sub getDependencyVersion {
    my $param = shift;

    if ( !exists $configuration->{$param} ) {
      die "Unknown parameter \"$param\"\n";
    }

    my $value = $configuration->{$param}->{'value'};
    my $version = "";
    my $tmpStr = "";
    if ( $param eq "RSCOUT_DIR" ) {
      # Simple case where the program reports it's version in it's usage details
      #./RepeatScout -blah
      #RepeatScout Version 1.0.6
      #
      #Usage: 
      #RepeatScout -sequence <seq> -output <out> -freq <freq> -l <l> [opts]
      #     -L # size of region to extend left or right (10000) 
      #     -match # reward for a match (+1)  
      #     -mismatch # penalty for a mismatch (-1) 
      #...
      $tmpStr = `$value/RepeatScout -blah 2>&1`;
      $tmpStr =~ /Version\s*(\d+\.\d+\.\d+)/;
      $version = $1;
    }elsif ( $param eq "RMBLAST_DIR" ) {
      #./rmblastn -version
      #rmblastn: 2.9.0+
      # Package: blast 2.9.0, build Sep  9 2019 15:21:42
      # The 3.x series prints a bare "rmblastn 3.0.7".
      $tmpStr = `$value/rmblastn -version`;
      $tmpStr =~ /rmblastn:?\s+(\d+\.\d+\.\d+\+?)/;
      $version = $1;
    }elsif ( $param eq "RECON_DIR" ) {
      # More complex.  The version is only printed in a
      # text file that may or may not have been installed by
      # the user. 
      # 00README:
      # Version 1.08 (Jan 2014)
      if ( -e "$value/00README" ) {
        open IN,"<$value/00README" or die;
        while (<IN>){
          if (/Version\s+(\d+\.\d+)/ )
          {
            $version = $1;
            last;
          }
        }
        close IN;
      }else {
        $tmpStr = `$value/eledef -v 2>&1`;
        $tmpStr =~ /RECON version\s+(\d\S+)/;
        $version = $1;
      }
    }elsif ( $param eq "TRF_DIR" ) { 
      # trf -V
      # Tandem Repeats Finder, Version 4.09
      # Copyright (C) Dr. Gary Benson 1999-2012. All rights reserved.
      $tmpStr = `$value/trf -V 2>&1`;
      $tmpStr =~ /Version\s+(\d+\.\d+)/;
      $version = $1;
    }elsif ( $param eq "UCSCTOOLS_DIR" ) { 
      # Currently Jim doesn't provide an option get a version.  Need to 
      # sneakily pull it out of the binaries strings.
      $tmpStr = `strings $value/twoBitInfo | grep "kent source version"`;
      $version = "unknown";
      $version = $1 if ( $tmpStr =~ /kent source version\s+(\d+)/ );
    }elsif ( $param eq "LTR_FORAGE_DIR" ) {
      # ltr-forage --version
      # ltr-forage 0.1.0
      $tmpStr = `$value/ltr-forage --version 2>&1`;
      $tmpStr =~ /ltr-forage\s+(\d+\.\d+\.\d+)/;
      $version = $1;
    }elsif ( $param eq "REPEATMASKER_DIR" ) {
      #RepeatMasker -v
      #RepeatMasker version 4.1.0
      $tmpStr = `$value/RepeatMasker -v 2>&1`;
      $tmpStr =~ /version\s+(\d+\.\d+\.\d+)/;
      $version = $1;
    }elsif ( $param eq "REPEATAFTERME_DIR" ) {
      #ExtendAlign -h
      #ExtendAlign Version 0.0.4 - build 436 date 20201118
      $tmpStr = `$value/RAMExtend -h 2>&1`;
      $tmpStr =~ /Version\s+(\d+\.\d+\.\d+)/;
      $version = $1;
    }elsif ( $param eq "SUFR_DIR" ) {
      #sufr --version
      #sufr 0.8.0
      $tmpStr = `$value/sufr --version 2>&1`;
      $tmpStr =~ /sufr\s+(\d+\.\d+\.\d+)/;
      $version = $1;
    }elsif ( $param eq "DFAM_CURATOR_DIR" ) {
      #te-composer --version
      #te-composer 0.6.6
      $tmpStr = `$value/te-composer --version 2>&1`;
      $tmpStr =~ /te-composer\s+(\d+\.\d+\.\d+)/;
      $version = $1;
    }
    return $version;
  }

}

1;
