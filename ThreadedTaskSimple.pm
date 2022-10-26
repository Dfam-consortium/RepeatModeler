#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) ThreadedTaskSimple.pm
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

ThreadedTaskSimple

=head1 SYNOPSIS

use ThreadedTaskSimple

Usage: 

    my $task = ThreadedTaskSimple->new();
    $task->addJob( name => "job-1", );

=head1 DESCRIPTION

A class for grouping a set of Job objects for a given multithreaded task.  
This is a simple implementation that blocks until execution is complete or an 
unrecoverable failure occurs.

=head1 INSTANCE METHODS

=cut

package ThreadedTaskSimple;
use strict;
use Cwd;
use Data::Dumper;
use Carp;
use ThreadedJob;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "ThreadedTaskSimple";

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

=head2 get_setDebug()

  Use: my $value    = getDebug();
  Use: my $oldValue = setDebug( $value );

  Get/Set the debug level.

=cut

##-------------------------------------------------------------------------##
sub getDebug
{
  my $obj = shift;

  my $value = $obj->{'debug'};

  return $value;
}

sub setDebug
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'debug'};
  $obj->{'debug'} = $value;

  return $oldValue;
}


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

=head2 get_setWorkingDir()

  Use: my $value    = getWorkingDir();
  Use: my $oldValue = setWorkingDir( $value );

  Get/Set the working dir.

=cut

##-------------------------------------------------------------------------##
sub getWorkingDir
{
  my $obj = shift;

  my $value = $obj->{'workingDir'};

  return $value;
}

sub setWorkingDir
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'workingDir'};
  $obj->{'workingDir'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setNumThreads()

  Use: my $value    = getNumThreads();
  Use: my $oldValue = setNumThreads( $value );

  Get/Set the number of threads.

=cut

##-------------------------------------------------------------------------##
sub getNumThreads
{
  my $obj = shift;

  my $value = $obj->{'numThreads'};

  return $value;
}

sub setNumThreads
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'numThreads'};
  $obj->{'numThreads'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setCompletionStatus()

  Use: my $value    = getCompletionStatus();
  Use: my $oldValue = setCompletionStatus( $value );

  Get/Set flag controlling the display (stdout) of completion
  status.

=cut

##-------------------------------------------------------------------------##
sub getCompletionStatus
{
  my $obj = shift;

  my $value = $obj->{'showCompletionStatus'};

  return $value;
}

sub setCompletionStatus
{
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'showCompletionStatus'};
  $obj->{'showCompletionStatus'} = $value;

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




sub addJob
{
  my $obj = shift;
  if ( @_ )
  {
    my $jb = ThreadedJob->new( @_ );
    push @{ $obj->{'jobs'} }, $jb;
  }
}

sub getEstimatedCompletionTime()
{
  my $obj = shift;
  
  my $avgRuntime = $obj->averageJobRuntime();
  my $remainingJobs = $obj->queuedJobsCount();
  my $runningJobs = $obj->runningJobsCount();
  if ( $avgRuntime <= 0 )
  {
    return ( -1 );
  }else {
    return( (($remainingJobs+$runningJobs) * $avgRuntime) / $obj->getNumThreads() );
  }
}

sub finishedJobsCount
{
  my $obj = shift;

  my $jobStatus = $obj->{'jobStatus'};

  my $jobsCnt = 0;
  foreach my $jobNum ( keys( %{$jobStatus} ) )
  {
    if ( $jobStatus->{$jobNum}->{'status'} eq "finished" )
    {
      $jobsCnt++;
    }
  }
  return $jobsCnt;
}

sub runningJobsCount
{
  my $obj = shift;

  my $jobStatus = $obj->{'jobStatus'};

  my $jobsCnt = 0;
  foreach my $jobNum ( keys( %{$jobStatus} ) )
  {
    if ( $jobStatus->{$jobNum}->{'status'} eq "running" )
    {
      $jobsCnt++;
    }
  }
  return $jobsCnt;
}


sub queuedJobsCount
{
  my $obj = shift;

  my $jobStatus = $obj->{'jobStatus'};

  my $jobsCnt = 0;
  foreach my $jobNum ( keys( %{$jobStatus} ) )
  {
    if ( $jobStatus->{$jobNum}->{'status'} eq "queued" )
    {
      $jobsCnt++;
    }
  }
  return $jobsCnt;
}

sub jobCount 
{
  my $obj = shift;
  return ( scalar(@{$obj->{'jobs'}}) ); 
}

sub averageJobRuntime
{
  my $obj = shift;

  my $jobStatus = $obj->{'jobStatus'};

  my $sumRuntimes = 0;
  my $finishedJobs = 0;
  foreach my $jobNum ( keys( %{$jobStatus} ) )
  {
    if ( $jobStatus->{$jobNum}->{'status'} eq "finished" )
    {
      my $runTime = $jobStatus->{$jobNum}->{'end_time'} - $jobStatus->{$jobNum}->{'start_time'};
      $sumRuntimes += $runTime;
      $finishedJobs++;
    }
  }
  if ( $finishedJobs == 0 )
  {
    return -1;
  }
  return ( int(sprintf("%d",$sumRuntimes / $finishedJobs)) );
}

sub percentComplete {
  my $obj = shift;
  return ( sprintf("%.1f %", ($obj->finishedJobsCount() / $obj->jobCount())*100 ));
}

sub isComplete
{
  my $obj = shift;

  my $jobStatus = $obj->{'jobStatus'};

  #print "" . Dumper($jobStatus) . "\n";

  my $flagComplete = 1;
  foreach my $jobNum ( keys( %{$jobStatus} ) )
  {
    if ( $jobStatus->{$jobNum}->{'status'} =~ /queued|running/ ) 
    {
      $flagComplete = 0;
      last;
    }
  }
  return $flagComplete;

}



sub printJobs {
  my $obj = shift;

  my $name = "unnamed_task";
  if ( defined $obj->getName()
       && $obj->getName() ne "" )
  {
    $name = $obj->getName();
  }
  my $jobStatus = $obj->{'jobStatus'};

  print "Task: $name, " . scalar(keys( %{$jobStatus} )) . " jobs\n";

  foreach my $jobNum ( sort { $a <=> $b } keys( %{$jobStatus} ) )
  {
    my $duration = $jobStatus->{$jobNum}->{'end_time'} - 
                   $jobStatus->{$jobNum}->{'start_time'};
    print "$jobNum: $jobStatus->{$jobNum}->{'name'} " . 
          "status=$jobStatus->{$jobNum}->{'status'} " .
          "retries=$jobStatus->{$jobNum}->{'retries'} "  .
          "retval=$jobStatus->{$jobNum}->{'retval'} " . 
          "start_time=".localtime($jobStatus->{$jobNum}->{'start_time'}) . 
          " end_time=" .localtime($jobStatus->{$jobNum}->{'end_time'}) .  
          " duration=$duration\n";
  }
}

##
## blocking_local_execute
##
sub execute
{
  my $obj = shift;

  my $workingDir = $obj->getWorkingDir();
  $workingDir = "." if ( ! defined $workingDir );

  my $maxJobs = $obj->getNumThreads();
  if ( !defined $maxJobs || $maxJobs < 1 )
  {
    croak "$CLASS"
        . "::local_execute(): NumThreads is not defined or is less than 1.\n";
  }
  my $maxRetries = $obj->getMaxRetries();
  $maxRetries = 0 if ( !defined $maxRetries );

  # Clear an initialize a hash to hold the runtime stats on the jobs
  $obj->{'jobStatus'} = {};
  my $jobNum    = 0;
  foreach my $job ( @{ $obj->{'jobs'} } )
  {
    $obj->{'jobStatus'}->{$jobNum}->{'name'}  = $job->getName();
    $obj->{'jobStatus'}->{$jobNum}->{'status'}  = 'queued';
    $obj->{'jobStatus'}->{$jobNum}->{'pid'}     = -1;
    $obj->{'jobStatus'}->{$jobNum}->{'retval'}  = undef;
    $obj->{'jobStatus'}->{$jobNum}->{'retries'} = 0;
    $obj->{'jobStatus'}->{$jobNum}->{'start_time'} = undef;
    $obj->{'jobStatus'}->{$jobNum}->{'end_time'} = undef;
    $jobNum++;
  }

  my %children     = ();
  my $badForkCount = 0;
  my $maxBadForks = 5; # TODO: parameterize
  my $failedJobs = 0;
  my $failEarly = 1;  # TODO: parameterize

  # Preserve original working directory
  my $origWorkingDir = getcwd(); 
  # Change to the preferred working directory for this task
  chdir($workingDir);

  # Start task manager loop
  while ( 1 )
  {

    # First take care of our children...like a good parent
    if ( keys( %children ) )
    {
      # Wait for at least one to exit;
      print "Waiting for a child to finish or die ( "
            . scalar( keys( %children ) )
            . " running )\n"
            if ( $obj->{'debug'} );
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
        $obj->{'jobStatus'}->{$jobNum}->{'pid'}    = -1;
        $obj->{'jobStatus'}->{$jobNum}->{'retval'} = $retVal; 

        # Check it's status
        if ( $retVal == 0 )
        {
          print "Child finished. PID=$childPID, "
              . "jobNum=$jobNum, RetVal=$retVal\n"
              if ( $obj->{'debug'} );
          $obj->{'jobStatus'}->{$jobNum}->{'status'} = 'finished';
          $obj->{'jobStatus'}->{$jobNum}->{'end_time'} = time();
          # TODO: Other validation?  Consider supporting a validation function
        } else
        {
          ## Process failed.
          print "Child die'd.  Organize the funeral for PID=$childPID, "
              . "jobNum=$jobNum, RetVal=$retVal\n"
              if ( $obj->{'debug'} );

          # Check how many times we have run it
          if ( $obj->{'jobStatus'}->{$jobNum}->{'retries'} < $maxRetries )
          {
            ## Under the retry limit...rerun it!
            print "Child for batch=$jobNum failed ($retVal) " . "retry#"
                . $obj->{'jobStatus'}->{$jobNum}->{'retries'} . "\n"
                if ( $obj->{'debug'} );
            $obj->{'jobStatus'}->{$jobNum}->{'retries'}++;
            $obj->{'jobStatus'}->{$jobNum}->{'status'} = 'queued';
            print "WARNING: Retrying job ( $jobNum ) [ $retVal ]...\n";
          } else
          {
            ## Too many retries.
            $obj->{'jobStatus'}->{$jobNum}->{'status'} = 'failed';
            $obj->{'jobStatus'}->{$jobNum}->{'end_time'} = localtime;
            $failedJobs++;
            if ( $failEarly ) {
              foreach my $jid ( keys(%{$obj->{'jobStatus'}}) ) {
                $obj->{'jobStatus'}->{$jid}->{'status'} = 'canceled' 
                     if ( $obj->{'jobStatus'}->{$jid}->{'status'} =~ /running|queued/ );
              }
              last; 
            }
          }
        }
      }
      # Report on completion if requested
      if ( exists $obj->{'showCompletionStatus'} && $obj->{'showCompletionStatus'} ) {
        my $cd = $obj->getEstimatedCompletionTime();
        if ( $cd > 0 ){
          my $minutes = sprintf( "%02s", int( $cd / 60 ) );
          $cd -= $minutes * 60;
          my $hours = sprintf( "%02s", int( $minutes / 60 ) );
          $minutes -= $hours * 60;
          my $seconds = sprintf( "%02s", int($cd) );
          print "     " . $obj->percentComplete() .  
                " complete. $hours:$minutes:$seconds (hh:mm:ss) est. time remaining\n";
        }
      }
    } # if ( keys %children...

    # Gather a list of batches to work on
    my @queuedJobs = grep { ( $obj->{'jobStatus'}->{$_}->{'status'} eq 'queued' ) }
        sort( { $a <=> $b } keys( %{$obj->{'jobStatus'}}) );

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
      print "Forking ( " . $job->getName() . " )\n" if ( $obj->{'debug'} );
    FORK:
      my $pid = fork();
      if ( $pid == 0 )
      {
        # Child
        print "Child Running " . $job->getName() . "\n" if ( $obj->{'debug'} );
        my $childFunction = $job->getFunction();
        my $childParameters = $job->getParameters();
        my $result = $childFunction->(@{$childParameters});
        exit( $result );
      } elsif ( $pid )
      {
        # Parent
        $children{$pid}                             = $queuedJobs[ $k ];
        $obj->{'jobStatus'}->{ $queuedJobs[ $k ] }->{'pid'}    = $pid;
        $obj->{'jobStatus'}->{ $queuedJobs[ $k ] }->{'status'} = 'running';
        $obj->{'jobStatus'}->{ $queuedJobs[ $k ] }->{'start_time'} = time();
      } elsif ( $! =~ /No more process/ )
      {
        $badForkCount++;
        last if ( $badForkCount > $maxBadForks );
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
  }
  # Restore working directory
  chdir($origWorkingDir);
  
  # Return stats ( TODO: Make this it's own function on the object data )
  my $cntQueued = 0;
  my $cntFailed = 0;
  my $cntRunning = 0;
  my $cntFinished = 0;
  foreach my $jobIdx ( keys(%{$obj->{'jobStatus'}}) )
  {
    my $status = $obj->{'jobStatus'}->{$jobIdx}->{'status'};
    if ( $status eq "queued" ) {
      $cntQueued++;
    }elsif ( $status eq "running" ) {
      $cntRunning++;
    }elsif ( $status eq "failed" ) {
      $cntFailed++;
    }elsif ( $status eq "finished" ) {
      $cntFinished++;
    }
  }
  my $taskFinished = 1;
  $taskFinished = 0 if ( $cntRunning > 0 || $cntFailed > 0 || $cntQueued > 0 );
    
  return( $taskFinished, $badForkCount, $cntFinished, $cntFailed );
}

1;
