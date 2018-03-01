#!/usr/bin/env perl
#
# mknoocnmap.pl                                    Erik Kluzek
#                                                  Dec/07/2011
#
# Create SCRIP grid and mapping files for a single-point or region
# that is assumed to be a land land-only region.
#
use Cwd;
use strict;
use English;
use IO::File;
use Getopt::Long;

#
# Global constants
#
my $degsiz = 0.1;

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


#-----------------------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]	    Gets map and grid files for a single land-only point.
REQUIRED OPTIONS
     -centerpoint [or -p] <lat,lon> Center latitude,longitude of the grid to create.
     -name [-or -n] <name>	    Name to use to describe point

OPTIONS
     -dx <number>                   Size of total grid in degrees in longitude direction 
                                    (default is $degsiz)
     -dy <number>                   Size of total grid in degrees in latitude direction 
                                    (default is $degsiz)
     -silent [or -s]		    Make output silent
     -help [or -h]		    Print usage to STDOUT.
     -verbose [or -v]		    Make output more verbose.
     -nx <number>                   Number of longitudes (default is 1)
     -ny <number>                   Number of latitudes  (default is 1)
EOF
}

#-----------------------------------------------------------------------------------------------

sub get_latlon {
#
# Return the latitude and longitude of the input string and validate it
#
  my $string = shift;
  my $desc   = shift;
  my $dx     = shift;
  my $dy     = shift;

  my $lat = undef;
  my $lon = undef;
  my $valreal1 = "[+-]?[0-9]*\.?[0-9]*[EedDqQ]?[0-9+-]*";

  if ( $string =~ /^($valreal1)\s*,\s*($valreal1)$/ ) {
     $lat = $1;
     $lon = $2;
  } else {
     die <<"EOF";
** $ProgName - Error in entering latitude/longitude for $desc **
EOF
  }
  if ( $dx <= 0.0 || $dx > 360. ) {
     die <<"EOF";
** $ProgName - Bad value for dx (=$dx) for $desc **
  }
  if ( $dy <= 0.0 || $dy > 180. ) {
     die <<"EOF";
** $ProgName - Bad value for dy (=$dy) for $desc **
  }
  if ( ($lat < -90.+$dy/2.0) || ($lat >  90.0-$dy/2.0) ) {
     die <<"EOF";
** $ProgName - Bad value for latitude (=$lat) for $desc **
EOF
  }
  if ( ($lon < $dx/2.0)   || ($lon > 360.0-$dx/2.0) ) {
     die <<"EOF";
** $ProgName - Bad value for longitude  (=$lat) for $desc **
EOF
  }
  return( $lat, $lon );

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

# Process command-line options

my %opts = (
            ctr     => undef,
            help    => undef,
            name    => undef,
            nx      => 1,
            ny      => 1,
            dx      => $degsiz,
            dy      => $degsiz,
            silent  => 0,
            verbose => 0,
           );

GetOptions(
    "p|centerpoint=s"  => \$opts{'ctr'},
    "n|name=s"         => \$opts{'name'},
    "nx=i"             => \$opts{'nx'},
    "ny=i"             => \$opts{'ny'},
    "dx=f"             => \$opts{'dx'},
    "dy=f"             => \$opts{'dy'},
    "h|help"           => \$opts{'help'},
    "s|silent"         => \$opts{'silent'},
    "v|verbose"        => \$opts{'verbose'},
)  or usage();

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

if ( $opts{'verbose'} && $opts{'silent'} ) {
    print "ERROR: Can NOT set both silent and verbose at once!\n";
    usage();
}
my $printlev;
if (      $opts{'verbose'} ) {
  $printlev = 2;
} elsif ( $opts{'silent'}  ) {
  $printlev = 0;
} else {
  $printlev = 1;
}

if ( ! defined($opts{'ctr'}) ) {
    print "ERROR: MUST set the center point\n";
    usage();
}
if ( ! defined($opts{'name'}) ) {
    print "ERROR: MUST set the name of the point\n";
    usage();
}
my $name = $opts{'name'};

my ($lat,$lon) = get_latlon( $opts{'ctr'}, $name, $opts{'dx'}, $opts{'dy'} );
my $S_lat = $lat - $opts{'dy'}/2.0;
my $N_lat = $lat + $opts{'dy'}/2.0;
my $W_lon = $lon - $opts{'dx'}/2.0;
my $E_lon = $lon + $opts{'dx'}/2.0;

my $nx = $opts{'nx'};
my $ny = $opts{'ny'};
if ( $opts{'nx'} < 1 ) {
    print "ERROR: nx MUST be greater than or equal to 1\n";
    usage();
}
if ( $opts{'ny'} < 1 ) {
    print "ERROR: ny MUST be greater than or equal to 1\n";
    usage();
}

#-----------------------------------------------------------------------------------------------
my $print;
if ( $printlev > 1 ) {
  $print = "PRINT=TRUE";  
}

# Creation date
my $cdate = `date +%y%m%d`; chomp( $cdate );

if ( $printlev > 0 ) {
  print "\n\nCreate SCRIP grid and mapping files for a single-point\n";
}
# land grid...
my $grddir  = absolute_path( "$scrdir/../mkmapgrids" );
my $grid1   = "$grddir/SCRIPgrid_${name}_nomask_c${cdate}.nc";
my $cmdenv  = "env S_LAT=$S_lat W_LON=$W_lon N_LAT=$N_lat E_LON=$E_lon " . 
             "NX=$nx NY=$ny PTNAME=$name $print ";

chdir( "$grddir" );
my $cmd     = "$cmdenv GRIDFILE=$grid1 ncl mkscripgrid.ncl";
if ( $printlev > 0 ) {
   print "Create land SCRIP gridfile\n";
   print "Execute: $cmd\n";
}
system( $cmd );

# ocean grid...
my $grid2 = "$grddir/SCRIPgrid_${name}_noocean_c${cdate}.nc";
my $cmd    = "$cmdenv GRIDFILE=$grid2 IMASK=0 ncl mkscripgrid.ncl";
if ( $printlev > 0 ) {
   print "Create ocean SCRIP gridfile\n";
   print "Execute: $cmd\n";
}
system( $cmd );

# Now create a unity mapping between the two...
# Note reversal of grid1 & grid2, because we want an ocean -> land
# mapping file
chdir( "$scrdir" );
my $mapfile = "map_${name}_noocean_to_${name}_nomask_aave_da_${cdate}.nc";
my $cmd = "env GRIDFILE1=$grid2 GRIDFILE2=$grid1 MAPFILE=$mapfile " .
          "$print ncl $scrdir/mkunitymap.ncl";

if ( $printlev > 0 ) {
  print "Create unity mapping file between the two gridfile\n";
  print "Execute: $cmd\n";
}
system( $cmd );

if ( $printlev > 0 ) {
  print "\n\nSuccessfully created grid/mapping files for single-point\n";
}
