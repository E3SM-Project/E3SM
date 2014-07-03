#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to return the decomposition information for CICE
#
# Usage:
#
# generate_cice_decomp [options]
#
# To get help on options and usage:
#
# generate_cice_decomp -help
#
#=======================================================================

use Cwd;
use strict;
#use diagnostics;
use Getopt::Long;
use English;

#-----------------------------------------------------------------------------------------------

#Figure out where configure directory is and where can use the XML/Lite module from
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives

my $cwd = getcwd();  # current working directory
my $cfgdir;
my $utilroot = $ENV{'UTILROOT'};

if ($ProgDir) { 
    $cfgdir = $ProgDir; 
} else { 
    $cfgdir = $cwd; 
}

# The XML::Lite module is required to parse the XML configuration files.
my $xmldir;
if (-f "../Tools/XML/Lite.pm") {
    $xmldir = "../Tools";
} elsif (-f "$cfgdir/../../../../scripts/ccsm_utils/Tools/perl5lib/XML/Lite.pm") {
    $xmldir = "$cfgdir/../../../../scripts/ccsm_utils/Tools/perl5lib";
} else {
    die <<"EOF";
** (generate_cice_decomp): Cannot find perl module \"XML/Lite.pm\" 
EOF
}
	  
#-----------------------------------------------------------------------------------------------
# Add $cfgdir/perl5lib to the list of paths that Perl searches for modules
my @dirs = ( $cfgdir, "$cfgdir/perl5lib", "$cfgdir/../../../../scripts/ccsm_utils/Tools/perl5lib", "$utilroot/Tools/perl5lib" );
unshift @INC, @dirs;
require XML::Lite;

my $result = eval "require Decomp::Config";
if ( ! defined($result) ) {
   die <<"EOF";
** (generate_cice_decomp): Cannot find perl module \"Decomp::Config\" from directories: @INC
EOF
}
require Decomp::Config;

#-----------------------------------------------------------------------------------------------
my $output = "all";
my $res = "gx1v6";
my $nx  = 320;
my $ny  = 384;

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]
OPTIONS
     -nproc <number>      (or -n)   Number of mpi tasks used.	
                                    (required)
     -res <resolution>    (or -r)   Horizontal resolution (gx1v6 etc.). 
                                    (default $res))
     -nx <number>                   number of lons 
                                    (required)
     -ny <number>                   number of lats
                                    (required)
     -thrds <number>      (or -t)   Number of threads per mpi task
                                    (default 1)
     -output <type>	  (or -o)   Either output: all, maxblocks, bsize-x, bsize-y, or decomptype
			            (default: $output)

EXAMPLES

   $ProgName -res gx1v6 -nx 320 -ny 384  -nproc 80 -output maxblocks

   will return a single value -- the optimum max number of blocks to use.

EOF
}

#------------------------------------------------------------------------------------------------

  my %opts = (
                res        => $res,
                nx         => $nx, 
                ny         => $ny, 
                nproc      => undef,
                thrds      => 1,
                output     => $output,
                printing   => 1,
                help       => 0,
                file       => "$cfgdir/cice_decomp.xml",
                spacecurve => 0,
                model      => "cice",
                platform   => "unset",
           );

  my $cmdline = @ARGV;
  GetOptions( 
              "r|res=s"      => \$opts{'res'},
              "nx=i"         => \$opts{'nx'},
              "ny=i"         => \$opts{'ny'},
              "n|nproc=i"    => \$opts{'nproc'},
              "t|thrds=i"    => \$opts{'thrds'},
              "o|output=s"   => \$opts{'output'},
              "h|elp"        => \$opts{'help'},
              "m|model=s"    => \$opts{'model'},     #no longer needed
              "s|spacecurve" => \$opts{'spacecurve'},#no longer needed
              "p|platform=s" => \$opts{'platform'},  #no longer needed
          ) or usage();

  # Check for unparsed arguments
  if (@ARGV) {
      print "ERROR: unrecognized arguments: @ARGV\n";
      usage();
  }
  if ( $opts{'help'} ) {
      usage();
  }

  foreach my $key ( keys( %opts ) ) {
     if ( $key ne "help" && ! defined($opts{$key}) ) {
        print "ERROR: required input $key was not set\n";
        usage();
     }
  }

$opts{'ProgName'} = $ProgName;
$opts{'ProgDir'}  = $cfgdir;
$opts{'cmdline'}  = $cmdline;

# Redefine nproc to be total procs, nproc*thrds
$opts{'nproc'} = $opts{'nproc'} * $opts{'thrds'};

# Set_horiz_grid sets the parameters for specific hgrid combinations.
my $nlat = $opts{'ny'};
my $nlon = $opts{'nx'};

# Try to read from the xml file
my $dcmp = Decomp::Config->new( \%opts );
my %decomp = ( maxblocks=>0, bsize_x=>0, bsize_y=>0, decomptype=>"", decompset=>"" );
my $matches = $dcmp->ReadXML( $opts{'file'}, \%decomp );

# If no xml entry, try to generate something
if ( $decomp{'maxblocks'} == 0) {
    %decomp = CalcDecompInfo( $nlat, $nlon, \%opts);
}

# adjust maxblocks to take into account threading
  $decomp{'maxblocks'} = $decomp{'maxblocks'} * $opts{'thrds'};

  if ( $decomp{'maxblocks'} == 0 ) {
     printf "%d %s",-1, "ERROR:($ProgName) No Decomp Created \n";
  } else {
     if (      $opts{'output'} eq "all"       ) {
	 printf "%d %d %d %d %d %s %s", 
	 $nlon, $nlat,
	 $decomp{'bsize_x'}, $decomp{'bsize_y'}, 
	 $decomp{'maxblocks'}, $decomp{'decomptype'}, $decomp{'decompset'};
      } elsif ( $opts{'output'} eq "maxblocks" ) {
        print $decomp{'maxblocks'};
      } elsif ( $opts{'output'} eq "bsize_x"   ) {
        print $decomp{'bsize_x'};
      } elsif ( $opts{'output'} eq "bsize_y"   ) {
        print $decomp{'bsize_y'};
      } elsif ( $opts{'output'} eq "decomptype") {
        print $decomp{'decomptype'};
      } elsif ( $opts{'output'} eq "decompset") {
        print $decomp{'decompset'};
      } else {
        print "ERROR:($ProgName) bad argument to output option $opts{'output'}\n";
        usage();
      }
      print "\n";
  }

#-----------------------------------------------------------------------------------------------
sub CalcDecompInfo {
#
# Calculate decomposition information
#  note that spacecurve must have nblocks in x and y direction divisible by only 2, 3, 5.
#  if spacecurve is the decomp target (dtypet) and a blocksize can't be found, use blkrobin
#  need to factor the nblocks to check if spacecurve blocksize is valid
#
  my $nlats    = shift;
  my $nlons    = shift;
  my $opts_ref = shift;

  my %opts   = %$opts_ref;
  my $nprocs = $opts{'nproc'};
  my $nthrds = $opts{'thrds'};

  my ($maxblocks,$bsize_x,$bsize_y,$decomptype,$decompset);
  my %decomp;
  my $set = 0;
  my $done = 0;
  my $nprocsx = 0;
  my $nprocsy = 0;
  my $nx = 0;
  my $ny = 0;
  my $nn = 0;
  my $nbx = 0;
  my $nby = 0;
  my $fac = 0;
  my $dtype;
  my $dtypet;
  my $mblocks = 0;
  my $bsize  = 0;
  my $bsizex = 0;
  my $bsizey = 0;
  my $nscore = 0.0 ;
  my $scok = 0;
  my $bscore = $nlons * $nlats * $nprocs ;

  $decomp{'decompset'} = "null";
  $decomp{'maxblocks'}  = 0;

  if ($nlats == 1) {
      $dtypet = "roundrobin";
      $dtype  = "roundrobin";
  } elsif ($nprocs * $nprocs > $nlons * $nlats * 6) {
#tcraig for testing  } elsif ($nprocs * $nprocs > 0) {
      $dtypet = "spacecurve";
      $dtype  = "blkrobin";
  } else {
      $dtypet = "blkrobin";
      $dtype  = "blkrobin";
  }

  if ($set == 0) {
     my $blksize = ($nlons * $nlats) / ($nprocs);
     $nn = 0;
     $done = 0;
     do {
         $nn = $nn + 1;
     } until ($nn * $nn > $blksize);

     my $nxnyfact = $nn * 4;
     if ($nlats == 1) {
	 $nxnyfact = $blksize;
     }

     #print "tcx1 $blksize $nn $nxnyfact \n";

     $blksize = ($nlons * $nlats) / ($nprocs * 8);
     $nx = 0;
     do {
         $nx = $nx + 1;
         if ($nlons % $nx == 0) {
             $ny = 0;
             do {
                 $ny = $ny + 1;
                 if ($nlats % $ny == 0) {
                     if ($nlats == 1) {
                        # min score is best score is max block size
                        $nscore = 1.0 / ($nx*$ny)
   		     } else {
                        # tcraig, somewhat arbitrary scoring system, best is min score
                        # min aspect ratio, match blksize, avoid very small blocks
                        $nscore = 0.5 * ($nx/$ny + $ny/$nx) + ($nx*$ny)/$blksize + $blksize/($nx*$ny) + 36/($nx*$ny);
                        $scok = 0;
                        if ($dtypet eq "spacecurve") {
                            # make sure nbx and nby are factorably by only 2, 3, 5
			    $nbx = $nlons / $nx;
                            $nby = $nlats / $ny;
                            $fac = 5;
                            #print "tcxf1 $nbx $nby\n";
                            do {
				if ($nbx % $fac == 0) {
				    $nbx = $nbx / $fac;
   			        } else {
 				    $fac = $fac - 1;
				}
                            } until ($fac == 1);
                            $fac = 5;
                            do {
				if ($nby % $fac == 0) {
				    $nby = $nby / $fac;
   			        } else {
 				    $fac = $fac - 1;
				}
                            } until ($fac == 1);
                            #print "tcxf2 $nbx $nby\n";
                            # if not, set nscore higher than bscore to skip it in search
                            if ($nbx == 1 && $nby == 1) {
                                $scok = 1;
			    }
                        }
		     }
                     #print "tcxc $nx $ny $nscore $bscore\n";
                     if ($dtypet eq "spacecurve" && $dtype ne "spacecurve" && $scok == 1) {
	                $bscore = $nscore;
                        $bsizex = $nx;
                        $bsizey = $ny;
                        $set    = 1;
                        $dtype = "spacecurve";
		     } elsif ($nscore < $bscore) {
                        if (($dtype eq "spacecurve" && $scok == 1) || ($dtype ne "spacecurve")) {
	                   $bscore = $nscore;
                           $bsizex = $nx;
                           $bsizey = $ny;
                           $set    = 1;
		        }
		     }
		 }
	     } until ($ny > $nxnyfact);
	 }
     } until ($nx >= $nxnyfact);
  }

  $bsize = $bsizex * $bsizey;

  #print "tcx0 $nlons $nlats $nprocs \n";
  #print "tcx1 $bsizex $bsizey $bscore \n";
  #print "tcx2 $bsize $blksize \n";

  if ($bsize > 0) {
      $mblocks = ($nlats * $nlons + $bsize * $nprocs - 1)/($bsize * $nprocs);
      $decomp{'nlats'}      = $nlats;
      $decomp{'nlons'}      = $nlons;
      $decomp{'bsize_x'}    = $bsizex;
      $decomp{'bsize_y'}    = $bsizey;
      $decomp{'maxblocks'}  = $mblocks;
      $decomp{'decompset'}  = "null";
      $decomp{'decomptype'} = $dtype;
  }

  return(%decomp);
}

#-------------------------------------------------------------------------------
sub clean
{
    my ($name) = @_;
    $name =~ s/^\s+//; # strip any leading whitespace 
    $name =~ s/\s+$//; # strip any trailing whitespace
    return ($name);
}

