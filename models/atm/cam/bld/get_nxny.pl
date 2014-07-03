#!/usr/bin/env perl

# get_nxny.sh reads nx and ny values from config_grid.xml in the CESM scripts.

use strict;
use Cwd;
use IO::File;
use IO::Handle;

my $grid_in = clean(@ARGV[0]);

# Find location of directory with this script.
(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
my $ProgDir = $1;                      # name of directory containing this script -- may be a
                                       # relative or absolute path, or null if the script is in
                                       # the user's PATH
my $cwd = getcwd();                    # current working directory
my $cfgdir;                            # absolute pathname of directory that contains this script
if ($ProgDir) { 
    $cfgdir = absolute_path($ProgDir);
} else {
    $cfgdir = $cwd;
}

# Add $cfgdir/perl5lib to the list of paths that Perl searches for modules
# and then get XML::Lite;
unshift @INC, "$cfgdir/perl5lib";
require XML::Lite;

# Find and open grid file
my $scriptsroot = "$cfgdir/../../../../scripts";
my $grid_file = "$scriptsroot/ccsm_utils/Case.template/config_grid.xml";
my $xml_grid = XML::Lite->new( $grid_file );

my @e = $xml_grid->elements_by_name( "gridhorz" );

my $nx, my $ny;

while ( my $e = shift @e ) {
    my %attr = $e->get_attributes();
    my $hgrid = clean($attr{'name'});
    my $alias = clean($attr{'alias'});
    if ($hgrid eq $grid_in or $alias eq $grid_in) {
        my @children = $e->get_children();
        foreach my $child (@children) {
            my $val  = clean($child->get_text());
            my $name = clean($child->get_name());
            if ($name eq 'nx') {$nx = $val;}
            if ($name eq 'ny') {$ny = $val;}
        }
        last;
    }
}

print "$nx $ny\n";

#-------------------------------------------------------------------------------
# Some subroutines copied from CESM scripts and configure.
#-------------------------------------------------------------------------------

sub clean
{
    my ($name) = @_;
    $name =~ s/^\s+//; # strip any leading whitespace 
    $name =~ s/\s+$//; # strip any trailing whitespace
    return ($name);
}

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
