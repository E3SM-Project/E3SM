#=======================================================================
#
# NAME
#
#   Decomp::OverallConfig -- A perl module to read in a CCSM4 Config XML file.
#
# SYNOPSIS
#
#  use Decomp::OverallConfig;
#
#  # Create a new decomp object
#  my $config = Decomp::OverallConfig->new( \%opts );
#  # Read in XML file with decomposition information for platform, res, #-procs, etc.
#  my %decomp = ();
#  my %decompresult = $config( $file, \%decomp );
#
# DESCRIPTION
#
#  This is a perl module to read in decomposition information from an XML file
#  for different: platforms, resolutions, #-processors, models etc.
#
# Methods:
#
#     new ---------------- Constructor.
#     ReadXML ------------ Read in the XML file and return a hash with data
#
# COLLABORATORS
# 
# XML::Lite
#
# HISTORY
#
# Date        Author                  Modification
#----------------------------------------------------------------------------------------
#  2007-Nov   Erik Kluzek             Original version
#----------------------------------------------------------------------------------------
#
#=======================================================================
use strict;
use XML::Lite;

package Decomp::OverallConfig;

#-------------------------------------------------------------------------------

sub new {
#
# Create a new DecompConfig object
#
# Required input to opts hash:
#
#     res ---------- resolution
#     mach --------- machine
#     compset ------ Compset
#     pecount ------ Approximate number of processors [S,M,L,XL]
#     ProgName ----- Main Program name
#     ProgDir ------ Main Program directory
#     cmdline ------ Main Program command line
#     printing ----- if debug printing should be done
#
  my $class     = shift;
  my $opts_ref  = shift;

  if ( ref($opts_ref) ne "HASH" ) { die "ERROR:: input opts is not a hash!\n"; }
  my %opts = %$opts_ref;

  # Error check that input and opts hash has the expected variables
  my $ProgName     = $opts{'ProgName'};
  my $nm = "${ProgName}::new";
  my @required_list = ( "res", "mach", "compset", "pecount", "cmdline",
                        "ProgDir", "printing" );
  my $self = {};
  foreach my $var ( @required_list ) {
     if ( ! defined($opts{$var}) ) {
        die "ERROR($nm): Required input variable $var was not found\n";
     }
     $self->{$var} = $opts{$var};
  }
  bless( $self, $class );
  return( $self );
}

#-------------------------------------------------------------------------------

sub ReadXML {
#
# Read in the XML file for the default decomposition configuration settings
#
# Input:
#
#   file --------- filename of XML file to read with a root of configInfo
#   decomp_ref --- Reference to a hash that has the expected output variables as keys
#
# Output:
#
#   decomp ------ Hash with same data as input decomp_ref with data values read from
#   tags in XML file that match (last match found).
#
# If variables are on the input XML file -- but NOT on the input decomp hash -- it
# will trigger an error and exit.
#
  my $self       = shift;
  my $file       = shift;
  my $decomp_ref = shift;

  my $ProgName = $self->{'ProgName'};
  my $nm = "${ProgName}::readXML";

  if ( ref($decomp_ref) ne "HASH" ) { die "ERROR::($nm) input decomp is not a hash!\n"; }

  # Initialize some local variables
  my $printing = $self->{'printing'};
  my $cmdline  = $self->{'cmdline'};
  my $matches  = undef;

  # Open file
#DBG  print "($nm) Read: $file\n" if $printing;
  my $xml = XML::Lite->new( $file );
  if ( ! defined($xml) ) {
    die "ERROR($nm): Trouble opening or reading $file\n";
  }
  #
  # Find the namelist element for this namelist
  #
  my $elm  = $xml->root_element( );
  my $root = "configInfo";
  my @list = $xml->elements_by_name( $root );
  if ( $#list < 0 ) {
    die "ERROR($nm): could not find the main $root root element in $file\n";
  }
  if ( $#list != 0 ) {
    die "ERROR($nm): $root root element in $file is duplicated, there should only be one\n";
  }
  $elm = $list[0];
  my @children = $elm->get_children();
  if ( $#children < 0 ) {
    die "ERROR($nm): There are no sub-elements to the $root element in $file\n";
  }
  #
  # Go through the sub-elements to the namelist element
  #
  $matches = 0;
  foreach my $child ( @children ) {
    #
    # Get the attributes for each element
    #
    my %atts = $child->get_attributes;
    # Name of element, and it's associated attributes
    my $name = $child->get_name();
    my @keys = keys(%atts);
    my $set = 1;
    if ( $#keys >= 0 ) {
      #
      # Check that all values match the appropriate settings
      #
      foreach my $key ( @keys ) {
#         printf "Key is %s\n", $key;
         foreach my $var ( ( "mach", "pecount" ) ) {
            # Match given var
            my $match = $atts{$key};
            if ( ($key eq $var) && ($self->{$var} !~ /^${match}$/ ) ) {
               $set = undef;
            }
         }
         # Match the begining of the string for the string for compset and res
         foreach my $var ( ( "compset", "res") ) {
            my $match = $atts{$key};
            if(($key eq $var) && ($self->{$var} !~ /^${match}/)  ) {
               $set = undef;
               # printf "match is %s\n", $match;
               # printf "self->{$var} is %s\n",$self->{$var};
            }
         }

      }
    }
    if ( $set ) {
       my @Grandchildren = $child->get_children();
       if ( $#Grandchildren > 0 ) {
          foreach my $Grandchild ( @Grandchildren ) {
             my $name  = $Grandchild->get_name();
             my $value = $Grandchild->get_text();
             if ( ! defined($$decomp_ref{$name}) ) {
                die "ERROR($nm): $name is NOT a valid element for the decomp output hash\n";
             }
             $$decomp_ref{$name} = $value;
          }
          $matches++;
       } else {
          die "ERROR($nm): No sub-elements to $name\n";
       }
    }
  }
#DBG  print "($nm): Matches = $matches\n" if $printing;
  return( $matches );
}

1 # To make use or require happy
