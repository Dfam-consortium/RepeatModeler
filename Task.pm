#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) Task.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A lightweight job management system for RepeatMasker/RepeatModeler.
##      This currently supports using the fork() command to run jobs locally
##      or the Sun Grid Engine job batching system for running on clusters.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2019 Developed by
#* Robert Hubley.
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

Task

=head1 SYNOPSIS

use Task

Usage: 

    my $jb = Task->new();

=head1 DESCRIPTION

A class for grouping a set of Job objects for a given task

=head1 INSTANCE METHODS

=cut

package Task;
use strict;
use Data::Dumper;
use Carp;
use Time::Piece;
use Job;
use Storable qw(nstore retrieve);
use LockFile::Simple qw(lock trylock unlock);

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "Task";
my $DEBUG = 1;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new
{
  my $class    = shift;
  my $arrayRef = shift;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

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
sub getName
{
  my $obj = shift;

  my $value = $obj->{'name'};

  return $value;
}

sub setName
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'name'};
  $obj->{'name'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setExecutor()

  Use: my $value    = getExecutor();
  Use: my $oldValue = setExecutor( $value );

  Get/Set the executor.

=cut

##-------------------------------------------------------------------------##
sub getExecutor
{
  my $obj = shift;

  my $value = $obj->{'executor'};

  return $value;
}

sub setExecutor
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'executor'};
  $obj->{'executor'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueue()

  Use: my $value    = getQueue();
  Use: my $oldValue = setQueue( $value );

  Get/Set the queue.

=cut

##-------------------------------------------------------------------------##
sub getQueue
{
  my $obj = shift;

  my $value = $obj->{'queue'};

  return $value;
}

sub setQueue
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'queue'};
  $obj->{'queue'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setMaxLocalJobs()

  Use: my $value    = getMaxLocalJobs();
  Use: my $oldValue = setMaxLocalJobs( $value );

  Get/Set the maxLocalJobs.

=cut

##-------------------------------------------------------------------------##
sub getMaxLocalJobs
{
  my $obj = shift;

  my $value = $obj->{'maxLocalJobs'};

  return $value;
}

sub setMaxLocalJobs
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'maxLocalJobs'};
  $obj->{'maxLocalJobs'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setMaxRetries()

  Use: my $value    = getMaxRetries();
  Use: my $oldValue = setMaxRetries( $value );

  Get/Set the maxRetries.

=cut

##-------------------------------------------------------------------------##
sub getMaxRetries
{
  my $obj = shift;

  my $value = $obj->{'maxRetries'};

  return $value;
}

sub setMaxRetries
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'maxRetries'};
  $obj->{'maxRetries'} = $value;

  return $oldValue;
}

sub job
{
  my $obj = shift;
  if ( @_ )
  {
    my $jb = Job->new( @_ );
    push @{ $obj->{'jobs'} }, $jb;
  }
}


##-------------------------------------------------------------------------##

=head2 execute()

  Use: $task->execute();

=cut

##-------------------------------------------------------------------------##
sub execute
{
  my $obj = shift;

  # Save state of jobs to disk
  # Decide where to execute
  if ( !defined $obj->{'executor'}
       || $obj->{'executor'} =~ /local/i )
  {
    $obj->_non_blocking_local_execute();
  } elsif ( $obj->{'executor'} =~ /sge/i )
  {
    if ( !defined $obj->{'queue'} )
    {
      croak "$CLASS"
          . "::execute(): Queue is not set.  A queue a required parameter for excutor SGE.\n";
    }
    $obj->_sge_execute();
  } else
  {
    croak "$CLASS" . "::execute(): Unknown executor \"$obj->{'executor'}\".\n";
  }

}


sub getEstimatedCompletionTime()
{
  my $obj = shift;
  
  my $avgRuntime = $obj->averageJobRuntime();
  my $remainingJobs = $obj->queuedJobsCount();
  if ( $avgRuntime <= 0 )
  {
    #
    return ( -1 );
  }else {
    return( $remainingJobs * $avgRuntime );
  }
}

sub queuedJobsCount()
{
  my $obj = shift;

  my $jobQueue = $obj->_queueStatus();

  my $jobsCnt = 0;
  foreach my $jobNum ( keys( %{$jobQueue} ) )
  {
    if ( $jobQueue->{$jobNum}->{'status'} eq "queued" )
    {
      $jobsCnt++;
    }
  }
  return $jobsCnt;
}

sub averageJobRuntime()
{
  my $obj = shift;

  my $jobQueue = $obj->_queueStatus();

  my $sumRuntimes = 0;
  my $finishedJobs = 0;
  foreach my $jobNum ( keys( %{$jobQueue} ) )
  {
    if ( $jobQueue->{$jobNum}->{'status'} eq "finished" )
    {
      my $runTime = $jobQueue->{$jobNum}->{'end_time'} - $jobQueue->{$jobNum}->{'start_time'};
      #print "Runtime: $runTime\n";
      $sumRuntimes += $runTime;
      $finishedJobs++;
    }
  }
  if ( $finishedJobs == 0 )
  {
    return -1;
  }
  return ( $sumRuntimes / $finishedJobs );
}

sub isComplete()
{
  my $obj = shift;

  my $jobQueue = $obj->_queueStatus();

  #print "" . Dumper($jobQueue) . "\n";

  my $flagComplete = 1;
  foreach my $jobNum ( keys( %{$jobQueue} ) )
  {
    if ( $jobQueue->{$jobNum}->{'status'} ne "finished" )
    {
      $flagComplete = 0;
      last;
    }
  }
  return $flagComplete;

}


## 
## Execute jobs using an Sun Grid Engine queue
##
sub _sge_execute
{
  my $obj = shift;

  my $date = localtime( time() );
  $date =~ s/[ ,\t,\n:]//g;
  my $uid  = "$$" . ".$date";
  my $name = "unnamed_task";
  if ( defined $obj->getName()
       && $obj->getName() ne "" )
  {
    $name = $obj->getName();
  }
  my $jobDataStore = "$name-$uid.dat";
  $obj->{'job_data_store'} = $jobDataStore;

  my %jobsQueue = ();
  my $jobNum = 0;
  foreach my $job ( @{ $obj->{'jobs'} } )
  {
    # TODO: Allow use of extra flags
    #       Set -o/-e  
    #       Set -S /bin/bash
    #       Turn command into script?
    my $qsubCmd = "qsub -q " . $obj->{'queue'} . " -N "
            . $job->getName()
            . " -b y -cwd "
            . $job->getCommand()
            . ""; 
    my $retVal = `$qsubCmd`;
    if ( $retVal =~ /Your job\s+(\d+)\s/ ) {
      $jobsQueue{$jobNum}->{'pid'}     = $1;
      $jobsQueue{$jobNum}->{'status'}  = 'queued';
    }else
    {
      # TODO: What do we do if we can't queue
    }
    $jobNum++;
  }
  &_saveDataStore( \%jobsQueue, $jobDataStore );
}

##
## _non_blocking_local_execute
##
sub _non_blocking_local_execute
{
  my $obj = shift;

  my $date = localtime( time() );
  $date =~ s/[ ,\t,\n:]//g;
  my $uid  = "$$" . ".$date";
  my $name = "unnamed_task";
  if ( defined $obj->getName()
       && $obj->getName() ne "" )
  {
    $name = $obj->getName();
  }
  my $jobDataStore = "$name-$uid.dat";
  $obj->{'job_data_store'} = $jobDataStore;

  my $maxJobs = $obj->getMaxLocalJobs();
  if ( !defined $maxJobs || $maxJobs < 1 )
  {
    croak "$CLASS"
        . "::local_execute(): MaxLocalJobs is not defined or is less than 1.\n";
  }
  my $maxRetries = $obj->getMaxRetries();
  $maxRetries = 0 if ( !defined $maxRetries );

  my %jobsQueue = ();
  my $jobNum    = 0;
  foreach my $job ( @{ $obj->{'jobs'} } )
  {
    $jobsQueue{$jobNum}->{'status'}  = 'queued';
    $jobsQueue{$jobNum}->{'pid'}     = -1;
    $jobsQueue{$jobNum}->{'retval'}  = undef;
    $jobsQueue{$jobNum}->{'retries'} = 0;
    $jobsQueue{$jobNum}->{'start_time'} = undef;
    $jobsQueue{$jobNum}->{'end_time'} = undef;
    $jobNum++;
  }

  &_saveDataStore( \%jobsQueue, $jobDataStore );

  # In order to be non-blocking we fork off a task manager that is
  # responsible for processing the job queue.  We then keep queue
  # state in a "Storable" file for asynchronous status checking.
  my $taskManagerPID = fork();
  if ( not $taskManagerPID )
  {

    #
    # TASK Manager Code
    #
    my %children     = ();
    my $badForkCount = 0;
    while ( 1 )
    {

      if ( keys( %children ) )
      {

        # Wait for at least one to exit;
        print "Waiting for a child to finish or die ( "
            . scalar( keys( %children ) )
            . " running )\n"
            if ( $DEBUG );
        my $childPID = wait();
        my $retVal   = ( $? >> 8 );

        # Check if we returned with a valid PID
        if ( $childPID > 0 )
        {
          ## Child process is gone
          # Find out what batch it was working on
          my $jobNum = $children{$childPID};

          # Delete it from the children list
          delete $children{$childPID};

          # Check it's status
          if ( $retVal == 0 )
          {
            print "Child finished. PID=$childPID, "
                . "jobNum=$jobNum, RetVal=$retVal\n"
                if ( $DEBUG );
            $jobsQueue{$jobNum}->{'pid'}    = -1;
            $jobsQueue{$jobNum}->{'status'} = 'finished';
            $jobsQueue{$jobNum}->{'end_time'} = localtime;
            # TODO: Other validation?
          } else
          {
            ## Process failed.
            print "Child die'd.  Organize the funeral for PID=$childPID, "
                . "jobNum=$jobNum, RetVal=$retVal\n"
                if ( $DEBUG );

            # Check how many times we have run it
            if ( $jobsQueue{$jobNum}->{'retries'} < $maxRetries )
            {
              ## Under the retry limit...rerun it!
              print "Child for batch=$jobNum failed ($retVal) " . "retry#"
                  . $jobsQueue{$jobNum}->{'retries'} . "\n"
                  if ( $DEBUG );
              $jobsQueue{$jobNum}->{'retries'}++;
              $jobsQueue{$jobNum}->{'pid'}    = -1;
              $jobsQueue{$jobNum}->{'status'} = 'queued';
              print "WARNING: Retrying job ( $jobNum ) [ $retVal ]...\n";
            } else
            {
              ## Too many retries.
              ## We are out of here!
              $jobsQueue{$jobNum}->{'pid'}    = -1;
              $jobsQueue{$jobNum}->{'status'} = 'failed';

              # TODO: Gracefully exit
            }
          }
          &_saveDataStore( \%jobsQueue, $jobDataStore );
        }
      }

      # Gather a list of batches to work on
      my @queuedJobs = grep { ( $jobsQueue{$_}->{'status'} eq 'queued' ) }
          sort( { $a <=> $b } keys( %jobsQueue ) );

      # Decide how many jobs to start
      my $numberToStart = scalar( @queuedJobs );
      $numberToStart = ( $maxJobs - keys( %children ) )
          if ( ( $maxJobs - keys( %children ) ) < $numberToStart );

      if ( @queuedJobs == 0 && !keys( %children ) )
      {
        # Looks like our work here is done
        last;
      }

      #
      # Loop through and fork to our hearts
      # content.
      #
      for ( my $k = 0 ; $k < $numberToStart ; $k++ )
      {
        my $job = $obj->{'jobs'}->[ $queuedJobs[ $k ] ];
        print "Forking ( " . $job->getName() . " )\n" if ( $DEBUG );
    FORK:
        my $pid = fork();
        if ( not $pid )
        {

          # Child
          print "Child Running " . $job->getName() . "\n";
          system( $job->getCommand() );
          exit;
        } elsif ( $pid )
        {

          # Parent
          $children{$pid}                             = $queuedJobs[ $k ];
          $jobsQueue{ $queuedJobs[ $k ] }->{'pid'}    = $pid;
          $jobsQueue{ $queuedJobs[ $k ] }->{'status'} = 'running';
          $jobsQueue{ $queuedJobs[ $k ] }->{'start_time'} = localtime;

          # TODO: Update job info on disk
        } elsif ( $! =~ /No more process/ )
        {
          $badForkCount++;
          sleep 5;
          redo FORK;
        } else
        {

          # Weird fork error
          print "\nWARNING: Cannot fork...for unknown reasons!\n";
          $badForkCount++;
          last;
        }
      }
      if ( $numberToStart )
      {
        &_saveDataStore( \%jobsQueue, $jobDataStore );
      }
    }
    exit;

  } elsif ( $taskManagerPID )
  {
    print "Task manager started (PID=$taskManagerPID). Executing jobs...\n";
  } elsif ( $! =~ /No more process/ )
  {

    # TODO: more sanity
    print "No more processes!!!!...\n";
  } else
  {

    # TODO: more sanity
    # Weird fork error
    print "\nWARNING: Cannot fork...for unknown reasons!\n";
  }
}

##
## _queueStatus
##
sub _queueStatus {
  my $obj = shift;

  my $jobQueue;
  if ( !defined $obj->{'executor'}
       || $obj->{'executor'} =~ /local/i )
  {
    $jobQueue = $obj->_localQueueStatus();
  } elsif ( $obj->{'executor'} =~ /sge/i )
  {
    $jobQueue = $obj->_SGEQueueStatus();
  } else
  {
    croak "$CLASS" . "::execute(): Unknown executor \"$obj->{'executor'}\".\n";
  }
  return $jobQueue;
}

##
## _localQueueStatus
##
sub _localQueueStatus {
  my $obj = shift;
  if ( !defined $obj->{'job_data_store'}
       || $obj->{'job_data_store'} eq "" )
  {
    return undef;
  }
  my $jobQueue = _retrieveDataStore( $obj->{'job_data_store'} );
  if ( scalar( keys( %{$jobQueue} ) ) != scalar( @{ $obj->{'jobs'} } ) )
  {
    # Either the count of jobs changed since the last execution or
    # something went wrong.  Either way it's "not complete".
    return undef;
  }

  return $jobQueue;
}

##
## _SGEQueueStatus
##
sub _SGEQueueStatus {
  my $obj = shift;
  
  if ( !defined $obj->{'job_data_store'}
       || $obj->{'job_data_store'} eq "" )
  {
    return undef;
  }
  my $jobQueue = _retrieveDataStore( $obj->{'job_data_store'} );
  if ( scalar( keys( %{$jobQueue} ) ) != scalar( @{ $obj->{'jobs'} } ) )
  {
    # Either the count of jobs changed since the last execution or
    # something went wrong.  Either way it's "not complete".
    return undef;
  }

  #
  # Update the job queue stats
  #
  # job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
  # -----------------------------------------------------------------------------------------------------------------
  # 50441 0.25000 rmblast-jo rhubley      r     03/14/2019 13:57:31 webserver@repeatmasker.systems     1        
  open QSTAT, "qstat |" or die "$CLASS" ."::isComplete(): Could not run qstat: $|\n";
  my %qstatStates = ();
  while ( <QSTAT> ) {
    if ( /^\s*(\d+)\s+\d+\.\d+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/ ) 
    {
      $qstatStates{$1} = $2;
      #my $timeStr = "$3 $4";
    }
  }
  close QSTAT;

  for ( my $i = 0; $i <= $#{$obj->{'jobs'}}; $i++ )
  {
    my $qjob = $jobQueue->{$i};
    my $stateCodes = $qstatStates{ $qjob->{'pid'} };
    if ( ! defined $stateCodes || $stateCodes eq "" ) 
    {
      # Check qacct
      my $qacctCmd = "qacct -j " . $qjob->{'pid'} . " 2>&1";
      my $status = `$qacctCmd`;
      if ( $status =~ /start_time\s+(\S+\s+\S+\s+\d+\s+\S+\s+\d\d\d\d)/ )
      { 
        # E.g Thu Mar 14 16:40:33 2019
        my $t = Time::Piece->strptime("$1", "%a %b %d %T %Y");
        $qjob->{'start_time'} = $t;
      }
      if ( $status =~ /end_time\s+(\S+\s+\S+\s+\d+\s+\S+\s+\d\d\d\d)/ )
      { 
        # E.g Thu Mar 14 16:40:33 2019
        my $t = Time::Piece->strptime("$1", "%a %b %d %T %Y");
        $qjob->{'end_time'} = $t;
      }
      if ( $status =~ /failed\s+(\d+)/ ) 
      {
        if ( $1 == 1 )
        {
          $qjob->{'status'} = 'failed';
          # TODO What else to do here
        }else {
          $qjob->{'status'} = 'finished';
        }
      }
      if ( $status =~ /exit_status\s+(\d+)/ ) 
      {
        $qjob->{'retval'} = $1;
      }
    }elsif ( $stateCodes =~ /r/ ) {
      $qjob->{'status'} = "running";
    }elsif ( $stateCodes =~ /q/ ) {
      $qjob->{'status'} = "queued";
    }
  }
  &_saveDataStore( $jobQueue, $obj->{'job_data_store'} );
  return $jobQueue;
}

##
## _retrieveDataStore
##
sub _retrieveDataStore
{
  my $filename = shift;

  my $limit = 100;
  for ( my $i = 0 ; $i < $limit ; $i++ )
  {
    if ( lock( $filename ) )
    {
      my $dataRef = retrieve( $filename );
      unlock( $filename );
      return $dataRef;
    }
  }
  return undef;
}

##
## _saveDataStore
##
sub _saveDataStore
{
  my $dataRef  = shift;
  my $filename = shift;

  my $limit = 100;
  for ( my $i = 0 ; $i < $limit ; $i++ )
  {
    if ( lock( $filename ) )
    {
      nstore( $dataRef, $filename );
      unlock( $filename );
      return 1;
    }
  }
  return 0;
}



##########################################
################OLD CODE##################
##########################################

sub blocking_local_execute
{
  my $obj = shift;

  my $jobDataStore = "task-name-date.dat";
  my $maxJobs      = $obj->getMaxLocalJobs();
  if ( !defined $maxJobs || $maxJobs < 1 )
  {
    croak "$CLASS"
        . "::local_execute(): MaxLocalJobs is not defined or is less than 1.\n";
  }
  my $maxRetries = $obj->getMaxRetries();
  $maxRetries = 0 if ( !defined $maxRetries );

  my %jobsQueue = ();
  my $jobNum    = 0;
  foreach my $job ( @{ $obj->{'jobs'} } )
  {
    $jobsQueue{$jobNum}->{'status'}  = 'queued';
    $jobsQueue{$jobNum}->{'pid'}     = -1;
    $jobsQueue{$jobNum}->{'retval'}  = undef;
    $jobsQueue{$jobNum}->{'retries'} = 0;
    $jobNum++;
  }
  &_saveDataStore( \%jobsQueue, $jobDataStore );

  my %children     = ();
  my $badForkCount = 0;
  while ( 1 )
  {

    if ( keys( %children ) )
    {

      # Wait for at least one to exit;
      print "Waiting for a child to finish or die ( "
          . scalar( keys( %children ) )
          . " running )\n"
          if ( $DEBUG );
      my $childPID = wait();
      my $retVal   = ( $? >> 8 );

      # Check if we returned with a valid PID
      if ( $childPID > 0 )
      {
        ## Child process is gone
        # Find out what batch it was working on
        my $jobNum = $children{$childPID};

        # Delete it from the children list
        delete $children{$childPID};

        # Check it's status
        if ( $retVal == 0 )
        {
          print "Child finished. PID=$childPID, "
              . "jobNum=$jobNum, RetVal=$retVal\n"
              if ( $DEBUG );
          $jobsQueue{$jobNum}->{'pid'}    = -1;
          $jobsQueue{$jobNum}->{'status'} = 'finished';

          # Other validation
        } else
        {
          ## Process failed.
          print "Child die'd.  Organize the funeral for PID=$childPID, "
              . "jobNum=$jobNum, RetVal=$retVal\n"
              if ( $DEBUG );

          # Check how many times we have run it
          if ( $jobsQueue{$jobNum}->{'retries'} < $maxRetries )
          {
            ## Under the retry limit...rerun it!
            print "Child for batch=$jobNum failed ($retVal) " . "retry#"
                . $jobsQueue{$jobNum}->{'retries'} . "\n"
                if ( $DEBUG );
            $jobsQueue{$jobNum}->{'retries'}++;
            $jobsQueue{$jobNum}->{'pid'}    = -1;
            $jobsQueue{$jobNum}->{'status'} = 'queued';
            print "WARNING: Retrying job ( $jobNum ) [ $retVal ]...\n";
          } else
          {
            ## Too many retries.
            ## We are out of here!
            $jobsQueue{$jobNum}->{'pid'}    = -1;
            $jobsQueue{$jobNum}->{'status'} = 'failed';

            # TODO: Gracefully exit
          }
        }
        &_saveDataStore( \%jobsQueue, $jobDataStore );
      }
    }

    # Gather a list of batches to work on
    my @queuedJobs = grep { ( $jobsQueue{$_}->{'status'} eq 'queued' ) }
        sort( { $a <=> $b } keys( %jobsQueue ) );

    # Decide how many jobs to start
    my $numberToStart = scalar( @queuedJobs );
    $numberToStart = ( $maxJobs - keys( %children ) )
        if ( ( $maxJobs - keys( %children ) ) < $numberToStart );

    if ( @queuedJobs == 0 && !keys( %children ) )
    {

      # Looks like our work here is done
      last;
    }

    #
    # Loop through and fork to our hearts
    # content.
    #
    for ( my $k = 0 ; $k < $numberToStart ; $k++ )
    {
      my $job = $obj->{'jobs'}->[ $queuedJobs[ $k ] ];
      print "Forking ( " . $job->getName() . " )\n" if ( $DEBUG );
  FORK:
      my $pid = fork();
      if ( not $pid )
      {

        # Child
        print "Child Running " . $job->getName() . "\n";
        system( $job->getCommand() );
        exit;
      } elsif ( $pid )
      {

        # Parent
        $children{$pid}                             = $queuedJobs[ $k ];
        $jobsQueue{ $queuedJobs[ $k ] }->{'pid'}    = $pid;
        $jobsQueue{ $queuedJobs[ $k ] }->{'status'} = 'running';

        # TODO: Update job info on disk
      } elsif ( $! =~ /No more process/ )
      {
        $badForkCount++;
        sleep 5;
        redo FORK;
      } else
      {

        # Weird fork error
        print "\nWARNING: Cannot fork...for unknown reasons!\n";
        $badForkCount++;
        last;
      }
    }
    if ( $numberToStart )
    {
      &_saveDataStore( \%jobsQueue, $jobDataStore );
    }
  }
  print "Done local_execute()\n";
}

1;
