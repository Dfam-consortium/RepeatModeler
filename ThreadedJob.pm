#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) ThreadedJob.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of a generic list collection utilizing
##      a perl array to hold data.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
###############################################################################
# ChangeLog
#
#     $Log$
#
###############################################################################
# To Do:
#
#

=head1 NAME

ThreadedJob

=head1 SYNOPSIS

use ThreadedJob

Usage: 

    my $jb = ThreadedJob->new();

=head1 DESCRIPTION

A class for storing a single job for a given task

=head1 INSTANCE METHODS

=cut

package ThreadedJob;
use strict;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "ThreadedJob";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class    = shift;
  my %nameValuePairs = @_;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  # Allow import of values
  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . uc( substr( $name, 0, 1 ) ) . substr( $name, 1 );
      eval { $this->$method( $value ); };
      if ( $@ ) {
        croak(
             "ThreadedTask::new Instance variable $name doesn't exist." . "" );
      }
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##


##-------------------------------------------------------------------------##

=head2 get_setName()

  Use: my $value    = getName();
  Use: my $oldValue = setName( $value );

  Get/Set the name.

=cut

##-------------------------------------------------------------------------##
sub getName {
  my $obj = shift;

  my $value = $obj->{'name'};

  return $value;
}

sub setName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'name'};
  $obj->{'name'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setFunction()

  Use: my $value    = getCommand();
  Use: my $oldValue = setCommand( $value );

  Get/Set the command.

=cut

##-------------------------------------------------------------------------##
sub getFunction {
  my $obj = shift;

  my $value = $obj->{'function'};

  return $value;
}

sub setFunction {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'function'};
  $obj->{'function'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setParameters()

  Use: my $values    = getParameters();
  Use: my $oldValues = setCpus( \@values );

  Get/Set the function parameters.

=cut

##-------------------------------------------------------------------------##
sub getParameters {
  my $obj = shift;

  my $value = $obj->{'parameters'};

  return $value;
}

sub setParameters {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'parameters'};
  $obj->{'parameters'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setWorkingDir()

  Use: my $value    = getWorkingDir();
  Use: my $oldValue = setWorkingDir( $value );

  Get/Set the workingDir.

=cut

##-------------------------------------------------------------------------##
sub getWorkingDir {
  my $obj = shift;

  my $value = $obj->{'workingDir'};

  return $value;
}

sub setWorkingDir {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'workingDir'};
  $obj->{'workingDir'} = $value;

  return $oldValue;
}

1;
