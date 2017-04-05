#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) dfamConsensusTool
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A command-line utitility for working with the Dfam_consensus
##      database.
##
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
# ChangeLog
#
#     $Log: dfamConsensusTool.pl,v $
#     Revision 1.4  2017/04/05 19:05:24  rhubley
#       - Doc improvements
#
#     Revision 1.3  2017/04/05 00:03:32  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

dfamConsensusTool - Command line tool for working with the Dfam_consensus database

=head1 SYNOPSIS

  dfamConsensusTool [options] -register   
  dfamConsensusTool [options] -validate <stockholm_file>
  dfamConsensusTool [options] -upload <stockholm_file>
  dfamConsensusTool [options] -classification

=head1 DESCRIPTION

The Dfam_consensus database is an open collection of Repetitive DNA consensus 
sequence models and corresponding seed alignments. It is directly compatible
with the RepeatMasker program and any consensus-based search tools. It is 
freely available and distribued under the Creative Commons Zero ( "CC0" ) 
license. The dfamConsensusTool.pl script is a command-line utility to aid 
with submission of new families to the Dfam_consensus database and is 
distributed as part of the RepeatModeler package. This utility provides the 
following basic features:

=over 4

Account Registration - Dfam_consensus submitters must have an 
account in order to submit families to the editors. The registration
process is quick and can be conducted at the website 
( http://www.repeatmasker.org/Dfam_consensus/#/login ) or through 
the dfamConsensus.pl tool itself.

Data Validation - The basic data format for Dfam/Dfam_conensus is 
a variant of the Stockholm format. This tool can be used to validate 
a particular Stockholm file prior to submission.

Data Submission - The primary role of this tool is to provide a 
reliable method for uploading a single curated family or a complete 
library of curated families to the database quickly and reliably.

=back

Full documentation can be found at: http://www.repeatmasker.org/RepeatModeler/dfamConsensusTool/

The options are:

=over 4

=item -version

Displays the version of the program

=item -register

Setup an account at the Dfam_consensus site.

=item -validate <stockholm_file>

Validate the format/data in a Stockholm file.

=item -upload  <stockholm_file>

Submit new families to the Dfam_consensus database.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2017 Robert Hubley, Institute for Systems Biology

=head1 LICENSE

This code may be used in accordance with the Open Source License v. 3.0
http://opensource.org/licenses/OSL-3.0

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Carp;
use Getopt::Long;
use Data::Dumper;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/..";
use MultAln;
use SeedAlignmentCollection;
use SeedAlignment;

## These are probably not already on the users
## default perl installation.
use JSON;
use URI;
use LWP::UserAgent;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name:  $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co Program" the
#  $Name:  $ macro will be blank so we should default
#  to what the ID tag for this file contains.
#
my $CVSNameTag = '$Name:  $';
my $CVSIdTag   =
    '$Id: dfamConsensusTool.pl,v 1.4 2017/04/05 19:05:24 rhubley Exp $';
my $Version = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG     = 0;
my $apiMajor  = 0;    # Version of the API we speak
my $apiMinor  = 0;    # ..
my $apiBugfix = 0;    # ..
my $dfamConsWebsite = "http://www.repeatmasker.org/Dfam_consensus";
my $server          = "www.repeatmasker.org";
my $port            = "10010";
my $adminEmail      = "help\@dfam.org";
my $failedFile;
my $username;
my $password;
my $sFilename;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',          # print out the version and exit
                    '-register',
                    '-upload=s',
                    '-validate=s',
                    '-classification',
                    '-status',
                    '-user=s',
                    '-pass=s',
                    '-server=s',
                    '-port=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) )
{
  usage();
}

sub usage
{
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} )
{
  print "$Version\n";
  exit;
}

$server = $options{'server'} if ( $options{'server'} );
$port   = $options{'port'}   if ( $options{'port'} );

unless (    $options{'validate'}
         || $options{'register'}
         || $options{'status'}
         || $options{'upload'}
         || $options{'classification'} )
{
  usage();
}

#
# Setup the User Agent
#
my $ua = LWP::UserAgent->new;
$ua->timeout( 10 );

#
# Make sure we have contact with the server and
# we can speak the same API version.
#
my $apiVersion = "";
my $res        = $ua->get( "http://$server:$port/version" );
if ( $res->decoded_content && !$res->is_error )
{
  my $data = from_json( $res->decoded_content );
  if ( !defined $data || !defined $data->{'major'} )
  {
    die "\nError in version response from the server [$server:$port]: "
        . Dumper( $res->decoded_content ) . "\n\n";
  }
  $apiVersion =
      $data->{'major'} . "." . $data->{'minor'} . "." . $data->{'bugfix'};

  if ( $data->{'major'} > $apiMajor )
  {
    die "\nError the server is running a newer API ( web: $apiVersion,\n"
        . "tool: $apiMajor.$apiMinor.$apiBugfix ) than this tool was designed for.\n"
        . "Please update RepeatModeler to obtain the most current version.\n\n";
  } elsif ( $data->{'minor'} > $apiMinor )
  {

    # TODO: Consider warning about minor level updates to API.
  }
} else
{
  print "\n[$server:$port - "
      . localtime()
      . "] Could not get a response from \n"
      . "the server.  Please try again, check the Dfam_consensus website\n"
      . "($dfamConsWebsite) for notices of scheduled maintenance or write\n"
      . "$adminEmail to report this problem.\n\n";
  die;
}

#
# Say hello
#
print "\nDfam_consensus Tool - version $Version ( API: $apiVersion )\n";
print "-------------------------------------------------------\n";

#
#  Register
#
if ( $options{'register'} )
{
  my $fullname;
  print "Registration to submit to Dfam_consensus is a two step process:\n"
      . "  1. Submit your account preferences using the form below.\n"
      . "  2. Send an email to $adminEmail to request \"submitter\" access\n"
      . "     to the system. Once this has been approved you will be able to\n"
      . "     begin uploading and viewing the status of sequences using this\n"
      . "     tool.\n";
  print "Email: ";
  chomp( $username = <STDIN> );
  print "Full Name: ";
  chomp( $fullname = <STDIN> );
  print "Password: ";
  system "stty -echo";
  chomp( $password = <STDIN> );
  system "stty echo";
  print "\n\n";

  if ( $username !~ /^(\S+)\@(\S+)$/ )
  {
    croak "Email address <$username> is not in the correct form.\n";
  }
  if ( $fullname eq "" )
  {
    croak "Missing Full Name.  This is a required field.\n";
  }
  if ( $password eq "" )
  {
    croak "Missing Password.  This is a required field.\n";
  }

  my $res = $ua->post(
                       "http://$server:$port/register",
                       {
                         'email'    => $username,
                         'name'     => $fullname,
                         'password' => $password
                       }
  );

  if ( $res->decoded_content || $res->is_error )
  {
    die "Error logging into server [$server:$port]: "
        . Dumper( $res->decoded_content ) . "\n";
  }

  print "Registration successful!  Now you can send an email\n"
      . "to \"$adminEmail\" and request \"submitter\" access to\n"
      . "Dfam_consensus.\n";
  exit;
}

#
# Validation part-1
#   - Parse the input file and attempt to do non-database
#     validation first. That way we do not bother the
#     user with login prompts just to do a pre-check of
#     the file.
#
my $stockholmFile = SeedAlignmentCollection->new();
if ( $options{'upload'} || $options{'validate'} )
{
  $sFilename = $options{'upload'};
  $sFilename = $options{'validate'} if ( $options{'validate'} );

  open my $IN, "<$sFilename"
      or die "Could not open up stockholm file $sFilename for reading!\n";
  $stockholmFile->read_stockholm( $IN );
  close $IN;
  print "File $sFilename contains " . $stockholmFile->size() . " families.\n";
  print "[test] Stockholm Format: valid\n";

}

#
#  Login
#
if ( $options{'validate'} )
{
  print "\nLogging into server [$server:$port] to further\n"
      . "validate the input file.  If you do not already have\n"
      . "an account you can register by running this tool with\n"
      . "the \"-register\" option.\n";
} else
{
  print "\nLogging into server [$server:$port]. If\n"
      . "you do not already have an account you can register\n"
      . "by running this tool with the \"-register\" option.\n";
}
if ( !$options{'user'} )
{
  print "Username: ";
  chomp( $username = <STDIN> );
} else
{
  $username = $options{'user'};
}

if ( !$options{'pass'} )
{
  print "Password: ";
  system "stty -echo";
  chomp( $password = <STDIN> );
  system "stty echo";
  print "\n\n";
} else
{
  $password = $options{'pass'};
}

#
# Obtain an OAuth token
#
$res = $ua->post(
                  "http://$server:$port/login",
                  {
                    'email'    => $username,
                    'password' => $password
                  }
);

my $data;
if ( $res->decoded_content && !$res->is_error )
{
  $data = from_json( $res->decoded_content );
  if ( !defined $data || !defined $data->{'token'} )
  {
    if ( defined $data && $data->{'message'} )
    {
      die "Error logging into server[$server:$port]: "
          . $data->{'message'} . "\n";
    } else
    {
      die "Error logging into server[$server:$port]: "
          . Dumper( $res->decoded_content ) . "\n";
    }
  }
} else
{
  print "[$server:$port - "
      . localtime()
      . "] Could not get a response from \n"
      . "the server.  Please try again, check the Dfam_consensus website\n"
      . "($dfamConsWebsite) for notices of scheduled maintenance or write\n"
      . "$adminEmail to report this problem.\n\n";
  die;
}

my $datestring = localtime();
print "Login at $datestring\n";

#
# Now place the token in the header of all future
# requests.
#
$ua->default_header( 'authorization' => "Bearer " . $data->{'token'} );

#
# Get the classification heirarchy
#
if ( $options{'classification'} )
{
  print "Dfam_consensus Classification System\n";
  print "  The following list represents the valid\n"
      . "  classes for families in Dfam_consensus:\n";
}
$res = $ua->get( "http://$server:$port/classification" );
my %validClassification = ();

if ( $res->decoded_content && !$res->is_error )
{

  #print ">" . $res->decoded_content . "<\n";
  my $data  = from_json( $res->decoded_content );
  my @stack = ( $data );
  while ( @stack )
  {
    my $node = pop( @stack );
    if ( !defined $node->{'found'} )
    {
      $node->{'found'} = 1;
      my $type = "";
      $type = "\"" . $node->{'repeatmasker_type'} . "\""
          if ( exists $node->{'repeatmasker_type'} );
      my $subtype = "";
      $subtype = ", \"" . $node->{'repeatmasker_subtype'} . "\""
          if ( exists $node->{'repeatmasker_subtype'} );

      #print "    " . $node->{'full_name'} . "\n"
      print "    \""
          . lc( $node->{'full_name'} )
          . "\" => [ $type $subtype ], \n"
          if ( $options{'classification'} );
      $validClassification{ lc( $node->{'full_name'} ) } = 1;
      foreach my $child ( @{ $node->{'children'} } )
      {
        push @stack, $child;
      }
    }
  }
} else
{
  print "[$server:$port - "
      . localtime()
      . "] Could not get a response from \n"
      . "the server.  Please try again, check the Dfam_consensus website\n"
      . "($dfamConsWebsite) for notices of scheduled maintenance or write\n"
      . "$adminEmail to report this problem.\n\n";

  #print "Error returned from LWP: " . $res->decoded_content . "\n";
  exit;
}
exit if ( $options{'classification'} );

#
# Validation revisited ( a.k.a part-2 )
#    - now do database-based validation and the slow consensus
#      generation.
#
if ( $options{'upload'} || $options{'validate'} )
{
  my $overallFailure = 0;

  # Gather some intel
  my %clades = ();
  my %ids    = ();
  for ( my $i = 0 ; $i < $stockholmFile->size() ; $i++ )
  {
    my $seedAlign = $stockholmFile->get( $i );
    my $id        = $seedAlign->getId();
    for ( my $j = 0 ; $j < $seedAlign->cladeCount() ; $j++ )
    {
      push @{ $clades{ $seedAlign->getClade( $j ) } }, $id;
    }
    $ids{$id}++;
  }

  # Sanity check the names
  my @repModelName = ();
  my @nonUniq      = ();
  my @tooLong      = ();
  foreach my $id ( keys %ids )
  {

    # rnd-1_family-138#LINE/L1
    if ( $id =~ /^rnd-\d+_family-\d+.*/ )
    {
      push @repModelName, $id;
    }
    if ( $ids{$id} > 1 )
    {
      push @nonUniq, $id;
    }
    if ( length( $id ) > 45 )
    {
      push @tooLong, $id;
    }
  }
  if ( @tooLong )
  {
    $overallFailure = 1;
    print "[test] Identifier Length: failed\n"
        . "  Identifiers are currently limited to 45 characters in length.\n"
        . "  The following identifiers are too long:\n";
    my $nameStr = join( ", ", @tooLong );
    $nameStr =~ s/(.{60}[^,]*,\s*)/$1\n/g;
    $nameStr = "      " . $nameStr;
    $nameStr =~ s/\n(\S.*)/\n      $1/g;
    $nameStr .= "\n" if ( substr( $nameStr, -1, 1 ) ne "\n" );
    print $nameStr;
  } else
  {
    print "[test] Identifier Length: valid\n";
  }
  if ( @repModelName )
  {
    $overallFailure = 1;
    print "[test] Custom Identifiers: failed\n"
        . "  Submissions to Dfam_consensus must be assigned names\n"
        . "  that are unique to the database.  Use of the RepeatModeler\n"
        . "  automatically assigned names would therefore not be\n"
        . "  unique.  Please see the help ( use \"-h\" ) for some\n"
        . "  pointers on naming your families and/or use the\n"
        . "  \"util/renameIDs.pl\" script to automatically rename\n"
        . "  the identifiers.\n\n"
        . "  The following identifiers look like a RepeatModeler\n"
        . "  assigned name:\n";
    my $nameStr = join( ", ", @repModelName );
    $nameStr =~ s/(.{60}[^,]*,\s*)/$1\n/g;
    $nameStr = "      " . $nameStr;
    $nameStr =~ s/\n(\S.*)/\n      $1/g;
    $nameStr .= "\n" if ( substr( $nameStr, -1, 1 ) ne "\n" );
    print $nameStr;
  } else
  {
    print "[test] Custom Identifiers: valid\n";
  }
  if ( @nonUniq )
  {
    $overallFailure = 1;
    print "[test] Unique (in-file) Identifiers: failed\n"
        . "  Non-unique seed alignment identifiers.  The following\n"
        . "  identifiers where used more than once in the stockholm\n"
        . "  file:\n";
    foreach my $id ( @nonUniq )
    {
      print "      $id\n";
    }
  } else
  {
    print "[test] Unique Identifiers: valid\n";
  }

  # Verify the clades
  my $failures = "";
  if ( keys( %clades ) )
  {
    foreach my $clade ( keys( %clades ) )
    {

      # Lookup in db
      my $url = URI->new( "http://$server:$port/taxonomy" );
      $url->query_form( { name => $clade, limit => 3 } );
      my $res = $ua->get( $url );
      if ( $res->decoded_content && !$res->is_error )
      {
        my $data = from_json( $res->decoded_content );
        if ( defined $data && defined $data->{'taxa'} )
        {
          if ( @{ $data->{'taxa'} } != 1
               || $data->{'taxa'}->[ 0 ]->{'species_name'} !~ /$clade/i )
          {
            $failures .=
                  "  Clade \"$clade\" found in seed alignments: "
                . join( ", ", @{ $clades{$clade} } ) . "\n"
                . "  does not uniquely map to an NCBI Taxonomy record.\n";
            if ( @{ $data->{'taxa'} } )
            {
              $failures .= "  Here are a few potential matches:\n";
              for my $name ( @{ $data->{'taxa'} } )
              {
                $failures .= "    " . $name->{'species_name'} . "\n";
              }
            }
          }
        } elsif ( defined $data && defined $data->{'message'} )
        {
          $failures .=
                "  Failed to verify clade \"$clade\" due to an\n"
              . "  error from the server:\n"
              . $data->{'message'} . "\n";
        } else
        {
          $failures .=
                "  Clade \"$clade\" found in seed alignments:\n" . "  "
              . join( ", ", @{ $clades{$clade} } ) . "\n"
              . "  does not map to an NCBI Taxonomy record.\n"
              . "  Please check the name at https://www.ncbi.nlm.nih.gov/taxonomy\n"
              . "  and correct the record in the stockholm file.\n";
        }
      } else
      {
        $failures .=
              "  Clade \"$clade\" found in seed alignments:\n" . "  "
            . join( ", ", @{ $clades{$clade} } ) . "\n"
            . "  does not map to an NCBI Taxonomy record.\n"
            . "  Please check the name at https://www.ncbi.nlm.nih.gov/taxonomy\n"
            . "  and correct the record in the stockholm file.\n";
      }
    }
    if ( $failures ne "" )
    {
      $overallFailure = 1;
      print "[test] NCBI Clade Names: failed\n$failures";
    } else
    {
      print "[test] NCBI Clade Names: valid\n";
    }
  }

# TODO: Verify classification
# TODO: Consider limiting the number of sequences per family -- don't overwhelm the db!
# TODO: Actually lookup IDs in db to make sure they are unique
# TODO: Lookup PMIDs?
# TODO: Verify DB Aliases
# TODO: Can we do a better job of validating the assembly/seqid/coordinates?

  if ( $overallFailure )
  {
    if ( $options{'validate'} )
    {
      print "\nFile failed validation.\n";
      exit;
    } else
    {
      print "\nFile failed validation - upload canceled.\n";
      exit;
    }
  } else
  {
    print "\nFile passes validation\n";
  }

  # Only used prefix of previously failed runs.  This way it will continue
  # to increment the ".failed.#" suffix.
  $failedFile = $sFilename;
  if ( $failedFile =~ /(.*)\.failed\.\d+$/ )
  {
    $failedFile = $1;
  }

  my $index = 1;
  while ( $index <= 15 && -s "$failedFile.failed.$index" )
  {
    $index++;
  }
  if ( $index == 16 )
  {
    die
"There are 15 version of the $failedFile.failed file in this directory.\n"
        . "Is there a problem we can help with?  Try emailing $adminEmail.\n";
  }
  $failedFile = "$failedFile.failed.$index";
}

if ( $options{'upload'} )
{

  # Basic request
  my $req = HTTP::Request->new( POST => "http://$server:$port/rmlibRepeat" );
  $req->content_type( 'application/json' );

  my $numErrs = 0;
  open OUT, ">$failedFile" or die "Could not open $failedFile for writing!\n";
  for ( my $i = 0 ; $i < $stockholmFile->size() ; $i++ )
  {
    my $seedAlign = $stockholmFile->get( $i );
    print "Working on " . $seedAlign->getId() . "\n";

    # Generate the consensus
    #   - TODO: Allow the users to use the RF line to specify their own???
    print "  - building consensus...\n";
    my @sequences = ();
    for ( my $j = 0 ; $j < $seedAlign->alignmentCount() ; $j++ )
    {
      my ( $assemblyName, $sequenceName, $start, $end, $orient, $sequence ) =
          $seedAlign->getAlignment( $j );

      # Replace prefix "." with spaces
      if ( $sequence =~ /^([\.]+)/ )
      {
        substr( $sequence, 0, length( $1 ) ) = " " x ( length( $1 ) );
      }

      # Replace suffixe "." with spaces
      if ( $sequence =~ /[^\.]([\.]+)$/ )
      {
        substr( $sequence, length( $sequence ) - length( $1 ), length( $1 ) ) =
            " " x ( length( $1 ) );
      }

      # Replace "." with "-"
      $sequence =~ s/\./\-/g;
      push @sequences, $sequence;
    }

    # Get the spaced (i.e with gaps) consensus
    my $consensus =
        MultAln::buildConsensusFromArray( sequences => \@sequences );
    $consensus =~ s/-//g;

    # Generate the minimal stockholm format to send along
    my $minStockholm = "# STOCKHOLM 1.0\n";
    $minStockholm .= "#=GC RF   " . $seedAlign->getRfLine() . "\n";
    for ( my $j = 0 ; $j < $seedAlign->alignmentCount() ; $j++ )
    {
      my ( $assemblyName, $sequenceName, $start, $end, $orient, $sequence ) =
          $seedAlign->getAlignment( $j );
      my $id;
      $id .= "$assemblyName:" if ( $assemblyName );
      $id .= "$sequenceName:" if ( $sequenceName );
      if ( $orient eq "+" )
      {
        $id .= "$start-$end";
      } else
      {
        $id .= "$end-$start";
      }
      $minStockholm .= "$id    $sequence\n";
    }
    $minStockholm .= "//";

    # prepare JSON
    my $record = {
                   'name'           => $seedAlign->getId(),
                   'sequence'       => $consensus,
                   'dfam_consensus' => "1",
                   'seed_alignment' => $minStockholm
    };

    # clades
    if ( $seedAlign->cladeCount() )
    {
      $record->{'clades'} = [];
      for ( my $j = 0 ; $j < $seedAlign->cladeCount() ; $j++ )
      {
        my $clade = $seedAlign->getClade( $j );
        push @{ $record->{'clades'} }, { 'name' => $clade };
      }
    }

    # description
    if ( $seedAlign->getDescription() )
    {
      $record->{'description'} = $seedAlign->getDescription();
    }

    # TODO: DB references

    # citations
    if ( $seedAlign->citationCount() )
    {
      $record->{'citations'} = [];
      for ( my $j = 0 ; $j < $seedAlign->citationCount() ; $j++ )
      {
        my ( $pmid, $title, $author, $journal ) = $seedAlign->getCitation( $j );
        my $citation = { 'pmid' => $pmid };
        $citation->{'title'}   = $title   if ( $title );
        $citation->{'author'}  = $author  if ( $author );
        $citation->{'journal'} = $journal if ( $journal );
        push @{ $record->{'citations'} }, $citation;
      }
    }

    # classification
    if ( $seedAlign->getClassification() )
    {
      $record->{'classification'} = $seedAlign->getClassification();
    }

    my $recJson = to_json( $record );

    print "  - uploading to the server...\n";

    #print " -- JSON = " . $recJson . "\n";
    # Send to the server
    $req->content( $recJson );
    my $res     = $ua->request( $req );
    my $message = "";
    if ( $res->decoded_content && !$res->is_error )
    {
      my $data = from_json( $res->decoded_content );

      #print "" . $res->decoded_content . "\n";
      if ( !defined $data || !defined $data->{'token'} )
      {
        $message = "[$server:$port - " . localtime() . "] Unknown failure";
        if ( defined $data && $data->{'message'} )
        {
          $message =
              "[$server:$port - " . localtime() . "] " . $data->{'message'};
          if ( defined $data->{'results'} )
          {
            $message .= " details: " . $data->{'results'};
          }
        }
      }
    }
    if ( $message ne "" )
    {
      print "   - upload failed: $message\n";

      # Write record + message to failed file.
      $seedAlign->setComment( $message . "\n" );
      print OUT "" . $seedAlign->toString();
      $numErrs++;
    }
  }
  close OUT;
  if ( $numErrs )
  {

    # TODO:
    if ( $numErrs == $stockholmFile->size() )
    {
      unlink( $failedFile )
          if ( -s $failedFile );
      print "\nNone of the families were uploaded successfully the server.\n"
          . "Please send the error messages to $adminEmail to seek help\n"
          . "with this matter.\n\n";
    } else
    {
      my $successful = $stockholmFile->size() - $numErrs;
      print "\n$successful out of "
          . $stockholmFile->size()
          . " families were uploaded\n"
          . "successfully.  Failed families have been copied to $failedFile\n"
          . "along with the error message for each.  Once these errors have been\n"
          . "corrected this file can be used to upload the failed families.\n\n";
    }
  }
} elsif ( $options{'status'} )
{
  print "Status is not implemented yet!\n";
}

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my _privateMethod( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##      Methods with the prefix "_" are conventionally considered
##      private.  This bit-o-documentation is not formatted to
##      print out when perldoc is run on this file.
##
##-------------------------------------------------------------------------##
sub _privateMethod
{
  my %parameters = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];
  print "$subroutine( " . @{ [ %parameters ] } . "): Called\n" if ( $DEBUG );

}

1;
