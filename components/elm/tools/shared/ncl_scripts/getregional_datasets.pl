#!/usr/bin/env perl
#=======================================================================
#
# Extract out regional datasets from the global datasets.
#
# Usage:
#
# getregional_datasets.pl
#
#  Erik Kluzek
#  Aug/28/2009
#  $Id: getregional_datasets.pl 55541 2013-11-22 08:41:13Z erik $
#  $HeadURL;
#
#=======================================================================

use Cwd;
use strict;
#use diagnostics;
use English;
use Getopt::Long;
use IO::File;

#-----------------------------------------------------------------------------------------------
# Set the directory that contains this scripts.  If the command was issued using a 
# relative or absolute path, that path is in $ProgDir.  Otherwise assume the
# command was issued from the current working directory.

(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
my $ProgDir = $1;                      # name of directory containing this script -- may be a
                                       # relative or absolute path, or null if the script is in
                                       # the user's PATH
my $cmdline = "@ARGV";                 # Command line arguments to script
my $cwd = getcwd();                    # current working directory
my $scrdir;                            # absolute pathname of directory that contains this script
my $nm = "ProgName::";                 # name to use if script dies
if ($ProgDir) { 
    $scrdir = absolute_path($ProgDir);
} else {
    $scrdir = $cwd;
}

my $gridfilename = "fatmlndfrc";

#-----------------------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]	    Extracts out files for a single box region from the global
                                    grid for the region of interest. Choose a box determined by
                                    the NorthEast and SouthWest corners.
REQUIRED OPTIONS
     -infilelist  "inlistfilename"  Input list of files to extract (namelist format).
       [-or -i] (REQUIRED)          The file variable ($gridfilename) is also required.
     -outfilelist "outlistfilename" Output list of filenames to extract (namelist format).
       [-or -o] (REQUIRED)          Must be the same variable names as infilelist.
     -NE_corner "lat,lon" [or -ne]  North East corner latitude and longitude                       (REQUIRED)
     -SW_corner "lat,lon" [or -sw]  South West corner latitude and longitude                       (REQUIRED)
OPTIONS
     -debug [or -d]                 Just debug by printing out what the script would do.
                                    This can be useful to find the size of the output area.
     -help [or -h]                  Print usage to STDOUT.

     -verbose [or -v]               Make output more verbose.
EOF
}

sub get_latlon {
#
# Return the latitude and longitude of the input string and validate it
#
  my $string = shift;
  my $desc   = shift;

  my $lat = undef;
  my $lon = undef;
  my $valreal1 = "[+-]?[0-9]*\.?[0-9]+[EedDqQ]?[0-9+-]*";

  if ( $string =~ /^($valreal1)\s*,\s*($valreal1)$/ ) {
     $lat = $1;
     $lon = $2;
  } else {
     die <<"EOF";
** $ProgName - Error in entering latitude/longitude for $desc **
EOF
  }
  if ( ($lat < -90.) || ($lat >  90.0) ) {
     die <<"EOF";
** $ProgName - Bad value for latitude (=$lat) for $desc **
EOF
  }
  if ( ($lon < 0.)   || ($lon > 360.0) ) {
     die <<"EOF";
** $ProgName - Bad value for longitude  (=$lat) for $desc **
EOF
  }
  return( $lat, $lon );

}

#-----------------------------------------------------------------------------------------------

# Process command-line options.

my %opts = ( 
              SW_corner        => undef,
              NE_corner        => undef,
              infilelist       => undef,
              outfilelist      => undef,
              help             => 0, 
              verbose          => 0,
              debug            => 0,
           );
GetOptions(
    "sw|SW_corner=s"   => \$opts{'SW_corner'},
    "ne|NE_corner=s"   => \$opts{'NE_corner'},
    "i|infilelist=s"   => \$opts{'infilelist'},
    "o|outfilelist=s"  => \$opts{'outfilelist'},
    "h|help"           => \$opts{'help'},
    "d|debug"          => \$opts{'debug'},
    "v|verbose"        => \$opts{'verbose'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

if ( ! defined($opts{'infilelist'}) || ! defined($opts{'outfilelist'}) ) {
    print "ERROR: MUST set both infilelist and outfilelist\n";
    usage();
}
if ( ! defined($opts{'SW_corner'}) || ! defined($opts{'NE_corner'}) ) {
    print "ERROR: MUST set both SW_corner and NE_corner\n";
    usage();
}

my ($S_lat,$W_lon) = get_latlon( $opts{'SW_corner'}, "SW" );
my ($N_lat,$E_lon) = get_latlon( $opts{'NE_corner'}, "NE" );

if ( $N_lat <= $S_lat ) {
    print "ERROR: NE corner latitude less than or equal to SW corner latitude\n";
    usage();
}
if ( $E_lon <= $W_lon ) {
    print "ERROR: NE corner longitude less than or equal to SW corner longitude\n";
    usage();
}

#-----------------------------------------------------------------------------------------------
my $debug;
if ( $opts{'debug'} ) {
  $debug = "DEBUG=TRUE";  
}
my $print;
if ( $opts{'verbose'} ) {
  $print = "PRINT=TRUE";  
}

my %infiles  = parse_filelist( $opts{'infilelist'}  );
my %outfiles = parse_filelist( $opts{'outfilelist'} );

(my $GRIDFILE, my $INFILES, my $OUTFILES) = get_filelists( \%infiles, \%outfiles );

write_usermods( \%outfiles );

my $cmd = "env S_LAT=$S_lat W_LON=$W_lon N_LAT=$N_lat E_LON=$E_lon " . 
          "GRIDFILE=$GRIDFILE OUTFILELIST=$OUTFILES INFILELIST=$INFILES " .
          "$debug $print ncl $scrdir/getregional_datasets.ncl";

print "Execute: $cmd\n";
system( $cmd );

#-------------------------------------------------------------------------------

sub parse_filelist {
#
# Parse a list of files (in "filename = 'filepath'" format) into a hash
#
  my $file = shift;

  # check that the file exists
  (-f $file)  or  die "$nm: failed to find filelist file $file";
  my $fh = IO::File->new($file, '<') or die "$nm: can't open file: $file\n";
  
  my %files = ( );
  my $valstring1 = '\'[^\']*\'';
  my $valstring2 = '"[^"]*"';
  while( my $line = <$fh> ) {
    if ( $line =~ m/^\s*(\S+)\s*=\s*($valstring1|$valstring2)$/ ) {
       my $var    = $1;
       my $string = $2;
       $string =~ s/'|"//g;
       if ( exists($files{$var}) ) {
          die "$nm: variable listed twice in file ($file): $var\n";
       }
       $files{$var} = $string;
    # Ignore empty lines or comments
    } elsif ( ($line =~ m/^\s*$/) || ($line =~ m/^\s*!/) ) {
       # ignore empty lines or comments
    } else {
       die "$nm: unexpected line in $file: $line\n";
    }
  }
  $fh->close; 

  return( %files );
}

#-------------------------------------------------------------------------------

sub get_filelists {
#
# Make sure file hashes compare correctly, and if so return in and out lists
#
  my $infiles_ref  = shift;
  my $outfiles_ref = shift;

  my @infiles  = sort( keys(%$infiles_ref ) );
  my @outfiles = sort( keys(%$outfiles_ref) );

  if ( $#infiles != $#outfiles ) {
     die "$nm: number of infiles is different from outfiles\n";
  }
  if ( "@infiles" ne "@outfiles" ) {
     die "$nm: list of infiles is different from outfiles list\n";
  }
  my $infilelist  = "";
  my $outfilelist = "";

  foreach my $file ( @infiles ) {
     my $infile  = $$infiles_ref{$file};
     if ( ! -f "$infile" ) {
        die "$nm: infile ($file) $infile does NOT exist!\n";
     }
     if ( $infilelist eq "" ) {
        $infilelist  = $infile;
     } else {
        $infilelist  .= ",$infile";
     }
     my $outfile = $$outfiles_ref{$file};
     if ( $outfilelist eq "" ) {
        $outfilelist  = $outfile;
     } else {
        $outfilelist  .= ",$outfile";
     }
     if ( -f "$outfile" ) {
        die "$nm: outfile ($file) $outfile already exists, delete it if you want to overwrite!\n";
     }
  }
  my $var = $gridfilename;
  my $gridfile = "";
  if ( exists($$infiles_ref{$var}) ) {
     $gridfile = $$infiles_ref{$var};
  } else {
      die "$nm: the grid file ($var) is required to be on the lists!\n";
  }

  return( $gridfile, $infilelist, $outfilelist );
}

#-------------------------------------------------------------------------------

sub write_usermods {
#
# Write the user_nl_clm and xmlchng_cmnds files out
# These can be used to setup a case after getregional_datasets is run.
#
  my $outfiles_ref = shift;

  my $cwd = getcwd();  # current working directory

  #
  # Write out the user_nl_clm file
  #
  my $usrnlfile = "user_nl_clm";
  my $fh = IO::File->new($usrnlfile, '>') or die "$nm: can't open file: $usrnlfile\n";

  my $outgridfile = undef;
  foreach my $file ( sort(keys(%$outfiles_ref)) ) {
     my $filepath = $$outfiles_ref{$file};
     # Add current directory on front of path if not an absolute path in filepath
     if ( $filepath !~ m/^\// ) {
        $filepath = "$cwd/$filepath";
     }
     # Write all filenames out besides the gridfilename
     if ( $file ne $gridfilename ) {
        print $fh "$file = '$filepath'\n";
     } else {
       $outgridfile = $filepath;
     }
  }
  $fh->close();
  #
  # Write out the xmlchnge_cmnds file
  #
  (my $filename = $outgridfile)=~ s!(.*)/!!;
  my $filedir   = $1;
  my $cmndsfile = "xmlchange_cmnds";
  my $fh = IO::File->new($cmndsfile, '>') or die "$nm: can't open file: $cmndsfile\n";
  print $fh "./xmlchange ATM_DOMAIN_PATH=$filedir\n";
  print $fh "./xmlchange LND_DOMAIN_PATH=$filedir\n";
  print $fh "./xmlchange ATM_DOMAIN_FILE=$filename\n";
  print $fh "./xmlchange LND_DOMAIN_FILE=$filename\n";
  $fh->close();
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

