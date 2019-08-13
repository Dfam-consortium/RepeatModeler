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
use Data::Dumper;
require Exporter;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );
@ISA         = qw(Exporter);
my $CLASS = "RepModelConfig";
BEGIN
{
##----------------------------------------------------------------------##
##     CONFIGURE THE FOLLOWING PARAMETERS FOR YOUR INSTALLATION         ##
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
          'ABBLAST_DIR' => {
                             'command_line_override' => 'abblast_dir',
                             'description' => 'The path to the installation of the ABBLAST
sequence alignment program.',
                             'environment_override' => 'ABBLAST_DIR',
                             'expected_binaries' => [
                                                      'blastp'
                                                    ],
                             'param_type' => 'directory',
                             'required' => 0,
                             'value' => '/usr/local/abblast'
                           },
          'CDHIT_DIR' => {
                           'command_line_override' => 'cdhit_dir',
                           'description' => 'The path to the installation of the CD-Hit
sequence clustering package.',
                           'environment_override' => 'CDHIT_DIR',
                           'expected_binaries' => [
                                                    'cd-hit'
                                                  ],
                           'param_type' => 'directory',
                           'required' => 0,
                           'value' => '/usr/local/cd-hit'
                         },
          'LTR_DETECTOR_DIR' => {
                                  'command_line_override' => 'ltr_detector_dir',
                                  'description' => 'The path to the installation of the LtrDetector
structural LTR repeatfinding program.',
                                  'environment_override' => 'LTR_DETECTOR_DIR',
                                  'expected_binaries' => [
                                                           'LtrDetector'
                                                         ],
                                  'param_type' => 'directory',
                                  'required' => 0,
                                  'value' => '/home/rhubley/src/LtrDetector/bin'
                                },
          'LTR_RETRIEVER_DIR' => {
                                   'command_line_override' => 'ltr_retriever_dir',
                                   'description' => 'The path to the installation of the LTR_Retriever
structural LTR analysis package.',
                                   'environment_override' => 'LTR_RETRIEVER_DIR',
                                   'expected_binaries' => [
                                                            'LTR_Retriever'
                                                          ],
                                   'param_type' => 'directory',
                                   'required' => 0,
                                   'value' => '/home/rhubley/src/LTR_retriever-2.6'
                                 },
          'MAFFT_DIR' => {
                           'command_line_override' => 'mafft_dir',
                           'description' => 'The path to the installation of the MAFFT
multiple alignment program.',
                           'environment_override' => 'MAFFT_DIR',
                           'expected_binaries' => [
                                                    'bin/mafft'
                                                  ],
                           'param_type' => 'directory',
                           'required' => 0,
                           'value' => '/usr/local/mafft'
                         },
          'NINJA_DIR' => {
                           'command_line_override' => 'ninja_dir',
                           'description' => 'The path to the installation of the Ninja
phylogenetic analysis package.',
                           'environment_override' => 'NINJA_DIR',
                           'expected_binaries' => [
                                                    'Ninja'
                                                  ],
                           'param_type' => 'directory',
                           'required' => 0,
                           'value' => '/home/rhubley/projects/NINJA/NINJA'
                         },
          'NSEG_PRGM' => {
                           'command_line_override' => 'nseg_prgm',
                           'description' => 'The full path including the name for the NSEG program.',
                           'environment_override' => 'NSEG_PRGM',
                           'param_type' => 'program',
                           'required' => 1,
                           'value' => '/usr/local/bin/nseg'
                         },
          'RECON_DIR' => {
                           'command_line_override' => 'recon_dir',
                           'description' => 'The path to the installation of the RECON
de-novo repeatfinding program.',
                           'environment_override' => 'RECON_DIR',
                           'expected_binaries' => [
                                                    'eledef',
                                                    'eleredef'
                                                  ],
                           'param_type' => 'directory',
                           'required' => 1,
                           'value' => '/usr/local/RECON'
                         },
          'REPEATMASKER_DIR' => {
                                  'command_line_override' => 'repeatmasker_dir',
                                  'description' => 'The path to the installation of RepeatMasker.',
                                  'environment_override' => 'REPEATMASKER_DIR',
                                  'expected_binaries' => [
                                                           'RepeatMasker'
                                                         ],
                                  'param_type' => 'directory',
                                  'required' => 1,
                                  'value' => '/usr/local/RepeatMasker'
                                },
          'RMBLAST_DIR' => {
                             'command_line_override' => 'rmblast_dir',
                             'description' => 'The path to the installation of the RMBLAST
sequence alignment program.',
                             'environment_override' => 'RMBLAST_DIR',
                             'expected_binaries' => [
                                                      'rmblastn'
                                                    ],
                             'param_type' => 'directory',
                             'required' => 1,
                             'value' => '/usr/local/rmblast-2.6.0/bin'
                           },
          'RSCOUT_DIR' => {
                            'command_line_override' => 'rscout_dir',
                            'description' => 'The path to the installation of the RepeatScout
de-novo repeatfinding program.',
                            'environment_override' => 'RSCOUT_DIR',
                            'expected_binaries' => [
                                                     'RepeatScout',
                                                     'build_lmer_table'
                                                   ],
                            'param_type' => 'directory',
                            'required' => 1,
                            'value' => '/usr/local/RepeatScout'
                          },
          'TRF_PRGM' => {
                          'command_line_override' => 'trf_prgm',
                          'description' => 'The full path including the name for the TRF program.',
                          'environment_override' => 'TRF_PRGM',
                          'param_type' => 'program',
                          'required' => 1,
                          'value' => '/usr/local/bin/trf409.linux64'
                        }
        };

  ## EDCFG --do-not-remove--

##                                                                      ##
##                      END CONFIGURATION AREA                          ##
##----------------------------------------------------------------------##
##----------------------------------------------------------------------##
##  Do not edit below this line                                         ##
##----------------------------------------------------------------------##

$DEBUGALL = 0;

# 
# Prompt for a specific parameter and update the object
#
sub promptForParam{
  my $param = shift;
  my $screenHdr = shift;

  if ( ! exists $configuration->{$param} ) {
    return;
  }

  # Grab defaults
  my $defaultValue;
  if ( $configuration->{$param}->{'param_type'} eq "directory" )
  {
    if ( exists $configuration->{$param}->{'expected_binaries'} &&
         @{$configuration->{$param}->{'expected_binaries'}} ) 
    {
      my $binary = $configuration->{$param}->{'expected_binaries'}->[0];
      $defaultValue = `/usr/bin/which $binary`;
      $defaultValue =~ s/[\n\r\s]+//g;
      $defaultValue =~ s/^(.*)\/$binary/$1/;
    }elsif ( exists $configuration->{$param}->{'value'} ) {
      $defaultValue = $configuration->{$param}->{'value'};
    }
  }elsif ( $configuration->{$param}->{'param_type'} eq "program" )
  {
    # The program type is used in cases where a single
    # script/binary is referenced and may not have the
    # exact name we expect.  TRF is a good example of this
    # as the binary is often distributed with names like:
    # trf409.linux64 etc..
    if ( exists $configuration->{$param}->{'value'} ) {
      $defaultValue = $configuration->{$param}->{'value'};
    }
  }

  my $value = "";
  my $validParam;
  do { 
    $validParam = 1;
    system("clear");
    if ( $screenHdr ) {
      print "$screenHdr\n";
    }else {
      print "\n\n\n\n";
    }
    print "" . $configuration->{$param}->{'description'} . "\n";

    # Prompt and get the value
    if ( $defaultValue ) {
      print "$param [$defaultValue]: ";
    }else {
      print "$param: ";
    }
    $value = <STDIN>;
    $value =~ s/[\n\r]+//g;
    if ( $value eq "" && $defaultValue )
    {
      $value = $defaultValue;
    }

    if ( $configuration->{$param}->{'param_type'} eq "directory" )
    {
      if ( -d $value ) {
        foreach my $binary ( @{$configuration->{$param}->{'expected_binaries'}} )
        {
          if ( ! -x "$value/$binary" )
          {
            print "\nCould not find the required program \"$binary\" inside\n" 
                . "the directory \"$value\"!\n\n";
            $validParam = 0;
            last;
          }elsif ( -d "$value/$binary" )
          {
            print "\nCould not find the required program \"$binary\" inside\n"
                . "the directory \"$value\"!  It appears to be the name of a\n"
                . "subdirectory.\n\n";
            $validParam = 0;
            last;
          }
        }
      }else { 
          print "\nCould not find the \"$value\" directory.\n\n";
          $validParam = 0;
      }   
    }elsif ( $configuration->{$param}->{'param_type'} eq "program" )
    {
      if ( ! -x $value ) {
        print "\nThe program \"$value\" doesn't appear to exist\n"
            . "or it's not executable!\n\n";
        $validParam = 0;
      }elsif ( -d $value ){
        print "\nThe value \"$value\" appears to be a directory rather\n" 
            . "than an executable binary or script!\n\n";
        $validParam = 0;
      }
    }

    if ( $validParam == 0 )
    {
      print "<PRESS ENTER TO CONTINUE, CTRL-C TO BREAK>\n";
      <STDIN>;
    }
  }while ( $validParam == 0 );
  $configuration->{$param}->{'value'} = $value;
}

#
# save values
#
sub updateConfigFile{
  open IN,"<$CLASS.pm" or die;
  open OUT,">new-$CLASS.pm" or die;
  my $inCfg;
  $Data::Dumper::Sortkeys = 1;
  while (<IN>){
    if ( /##\s+STCFG/ ) {
      $inCfg = 1;
      print OUT;
      my $cStr = Dumper($configuration);
      $cStr =~ s/\$VAR1/  \$configuration/;
      print OUT "$cStr\n";
    }elsif ( /##\s+EDCFG/ ) {
      $inCfg = 0;
      print OUT;
    }elsif ( ! $inCfg ) {
      print OUT;
    }
  }
  close IN;
  close OUT;
  rename("$CLASS.pm", "$CLASS.pm.bak");
  rename("new-$CLASS.pm", "$CLASS.pm");
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
      push @options, "-" . $configuration->{$param}->{'command_line_override'} . "=s";
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
      $pod_str .= "=item -" . $configuration->{$param}->{'command_line_override'} . " <string>\n\n";
      $pod_str .= $configuration->{$param}->{'description'} . "\n\n";
    }
  }
  if ( $pod_str ) {
    return( "\n=over 4\n\n" . $pod_str . "=back\n\n" );
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
    if ( exists $options->{$configuration->{$param}->{'command_line_override'}} )
    {
      $configuration->{$param}->{'value'} = $options->{$configuration->{$param}->{'command_line_override'}};
    } elsif ( exists $ENV{$configuration->{$param}->{'environment_override'}} )
    {
      $configuration->{$param}->{'value'} = $ENV{$configuration->{$param}->{'environment_override'}};
    }
  }
}


}

1;