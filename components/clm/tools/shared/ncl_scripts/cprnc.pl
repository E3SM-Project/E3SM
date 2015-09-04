#!/usr/bin/env perl
#=======================================================================
#
# Compare two NetCDF files that do NOT have a time axis.
#
# Usage:
#
# cprnc_nontimefile.pl file1 file2
#
#  Erik Kluzek
#  May/04/2011
#  $Id: getregional_datasets.pl 25177 2010-10-16 05:12:30Z erik $
#  $HeadURL;
#
#=======================================================================

use Cwd;
use strict;
use English;
use Getopt::Long;
use IO::File;

#-----------------------------------------------------------------------------------------------
# Set the directory that contains this scripts.  If the command was issued using a 
# relative or absolute path, that path is in $ProgDir.  Otherwise assume the
# command was issued from the current working directory.

(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
my $ProgDir = $1;                      # name of directory containing this script -- may be a
                                       # relative or absolute path, or null if the script
                                       # is in
                                       # the user's PATH
my $cmdline = "@ARGV";                 # Command line arguments to script
my $cwd = getcwd();                    # current working directory
my $scrdir;                            # absolute pathname of directory that contains this script
my $nm = "${ProgName}::";              # name to use if script dies
if ($ProgDir) { 
    $scrdir = absolute_path($ProgDir);
} else {
    $scrdir = $cwd;
}

#--------------------------------------------------------------------------

sub absolute_path {
#
# Convert a pathname into an absolute pathname, expanding any . or .. characters.
# Assumes pathnames refer to a local filesystem.
# Assumes the directory separator is "/".
#
  my $path = shift;
  my $cwd = getcwd();  # current working directory
  my $abspath;         # resulting absolute pathname

# Strip off any leading or trailing whitespace.  (This pattern won't match if
# there's embedded whitespace.
  $path =~ s!^\s*(\S*)\s*$!$1!;

# Convert relative to absolute path.

  if ($path =~ m!^\.$!) {          # path is "."
      return $cwd;
  } elsif ($path =~ m!^\./!) {     # path starts with "./"
      $path =~ s!^\.!$cwd!;
  } elsif ($path =~ m!^\.\.$!) {   # path is ".."
      $path = "$cwd/..";
  } elsif ($path =~ m!^\.\./!) {   # path starts with "../"
      $path = "$cwd/$path";
  } elsif ($path =~ m!^[^/]!) {    # path starts with non-slash character
      $path = "$cwd/$path";
  }

  my ($dir, @dirs2);
  my @dirs = split "/", $path, -1;   # The -1 prevents split from stripping trailing nulls
                                     # This enables correct processing of the input "/".

  # Remove any "" that are not leading.
  for (my $i=0; $i<=$#dirs; ++$i) {
      if ($i == 0 or $dirs[$i] ne "") {
          push @dirs2, $dirs[$i];
      }
  }
  @dirs = ();

  # Remove any "."
  foreach $dir (@dirs2) {
      unless ($dir eq ".") {
          push @dirs, $dir;
      }
  }
  @dirs2 = ();

  # Remove the "subdir/.." parts.
  foreach $dir (@dirs) {
    if ( $dir !~ /^\.\.$/ ) {
        push @dirs2, $dir;
    } else {
        pop @dirs2;   # remove previous dir when current dir is ..
    }
  }
  if ($#dirs2 == 0 and $dirs2[0] eq "") { return "/"; }
  $abspath = join '/', @dirs2;
  return( $abspath );
}

#--------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options] file1.nc file2.nc   Compares two NetCDF files
                                             and reports if the contents
                                             are identical or different.

OPTIONS

     -help [or -h]                 Print usage to STDOUT.
     -[no]break                    Switch on [off] break on first difference
                                   (default is -break, abort after first difference found)
     -verbose [or -v]              Make output more verbose.

EOF
}

# Process command-line options.

my %opts = ( 
              help             => 0, 
              break            => 1,
              verbose          => 0,
           );
GetOptions(
    "h|help"           => \$opts{'help'},
    "break!"           => \$opts{'break'},
    "v|verbose"        => \$opts{'verbose'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check that two filenames are given
if ($#ARGV != 1) {
    print "ERROR: Must provide two filenames\n";
    usage();
}
my $file1 = shift(@ARGV);
my $file2 = shift(@ARGV);

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

print "Compare $file1 to $file2\n";

if ( ! -f $file1 ) {
   print "ERROR: $file1 does NOT exist\n";
   usage();
}
if ( ! -f $file2 ) {
   print "ERROR: $file2 does NOT exist\n";
   usage();
}

#-----------------------------------------------------------------------------------------------
my $print;
if ( $opts{'verbose'} ) {
  $print = "PRINT=TRUE";  
}

my $break;
if ( $opts{'break'} ) {
   $break = "BREAKONDIFF=TRUE";
}

my $cmd = "env MYFILE1=$file1 MYFILE2=$file2 $print $break ncl $scrdir/cprnc.ncl";

print "Execute: $cmd\n";
system( $cmd );
