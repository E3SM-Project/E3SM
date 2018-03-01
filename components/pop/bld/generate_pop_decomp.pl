#!/usr/bin/env perl

#=======================================================================
#  This is a script to return the decomposition information for POP
#
# Usage:
#      generate_pop_decomp [options]
# To get help on options and usage:
#      generate_pop_decomp -help
#=======================================================================

use strict;
use Getopt::Long;
use English;
use POSIX;

#-----------------------------------------------------------------------------------------------
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program

#-----------------------------------------------------------------------------------------------
sub usage {
    die <<EOF;
SYNOPSIS
    $ProgName [options]

OPTIONS
    -ccsmroot <path>               Full pathname for ccsmroot 
                                   (optional, default is .)
    -nproc <number>      (or -n)   Number of processors to use
                                   (required)
    -res <resolution>    (or -r)   Horizontal resolution             
                                   (optional, default is gx1v6)
    -thrds <number>      (or -t)   Number of threads per processor   
                                   (optional, default 1)
    -output <type>	 (or -o)   Either output: all, maxblocks, bsize-x, bsize-y, or decomptype
                                   (optional: default is all)

EXAMPLES

   $ProgName -res gx1v6 -nproc 80 -output maxblocks

   will return a single value -- the optimum max number of blocks to use

EOF
}

#------------------------------------------------------------------------------------------------
my %opts = (
    ccsmroot => undef,
    nproc    => undef,
    res      => 'gx1v6',
    thrds    => 1,
    output   => 'all',
    printing => 1,
    help     => 0,
    );

GetOptions( 
    "ccsmroot=s"   => \$opts{'ccsmroot'},
    "r|res=s"      => \$opts{'res'},
    "n|nproc=i"    => \$opts{'nproc'},
    "t|thrds=i"    => \$opts{'thrds'},
    "o|output=s"   => \$opts{'output'},
    "h|elp"        => \$opts{'help'},
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

#------------------------------------------------------------------------------------------------
# Add perl5lib to the list of paths that Perl searches for modules
my $cesmroot = $opts{'ccsmroot'};
my @dirs = ("$cesmroot/cime/utils/perl5lib");
unshift @INC, @dirs;
require Decomp::Config;

#------------------------------------------------------------------------------------------------
# redefine nproc to be total procs, nproc*thrds
$opts{'nproc'} = $opts{'nproc'} * $opts{'thrds'};

# try to read from the xml file
my $dcmp = Decomp::Config->new( \%opts );

my %decomp = ( maxblocks=>0, bsize_x=>0, bsize_y=>0, decomptype=>"",
               nlats=>0, nlons=>0, nx_blocks=>0, ny_blocks=>0 );

my $file = "$cesmroot/components/pop/bld/pop_decomp.xml";
my $matches = $dcmp->ReadXML( $file, \%decomp );

# if no xml entry, try to generate something
if ( $decomp{'maxblocks'} == 0) {
   %decomp = CalcDecompInfo( $decomp{'nlats'}, $decomp{'nlons'}, \%opts);
}
 
# adjust maxblocks to take into account threading
$decomp{'maxblocks'} = $decomp{'maxblocks'} * $opts{'thrds'};

if ( $decomp{'maxblocks'} == 0 ) {
   printf "%d %s",-1, "ERROR:($ProgName) No Decomp Created \n";
} else {
   if (      $opts{'output'} eq "all"       ) {
     printf "%d %d %d %d %d %s %d %d", $decomp{'nlons'}, $decomp{'nlats'}, 
        $decomp{'bsize_x'}, $decomp{'bsize_y'}, $decomp{'maxblocks'}, $decomp{'decomptype'},
        $decomp{'nx_blocks'}, $decomp{'ny_blocks'};
    } elsif ( $opts{'output'} eq "maxblocks" ) {
      print $decomp{'maxblocks'};
    } elsif ( $opts{'output'} eq "bsize_x"   ) {
      print $decomp{'bsize_x'};
    } elsif ( $opts{'output'} eq "bsize_y"   ) {
      print $decomp{'bsize_y'};
    } elsif ( $opts{'output'} eq "decomptype") {
      print $decomp{'decomptype'};
    } elsif ( $opts{'output'} eq "nx_blocks") {
      print $decomp{'nx_blocks'};
    } elsif ( $opts{'output'} eq "ny_blocks") {
      print $decomp{'ny_blocks'};
    } else {
      print "ERROR:($ProgName) bad argument to output option $opts{'output'}\n";
      usage();
    }
    print "\n";
}

exit(0);

#------------------------------------------------------------------------------------------------
sub CalcDecompInfo {
    #
    # Calculate decomposition information
    # Tries to first find an even cartesian decomposition (set = 1)
    # If can't find an even decomp, tries a space filling curve (set = 2)
    # If it can't find spacecurve, find a blockone decomp
    #
    my $nlats    = shift;
    my $nlons    = shift;
    my $opts_ref = shift;
    
    my %opts   = %$opts_ref;
    my $nprocs = $opts{'nproc'};
    
    #  my ($maxblocks,$bsize_x,$bsize_y,$decomptype);
    my %decomp;
    my $found = 0;
    my $done = 0;
    my $nx = 0;
    my $ny = 0;
    my $nxblocks = 0;
    my $nyblocks = 0;
    my $nn = 0;
    my $tmp = 0;
    my $nscore = 0.0 ;
    my $bscore = $nlons * $nlats * $nprocs ;
    my $bsizex = 0;
    my $bsizey = 0;
    my $dtype;
    my $mblocks = 0;

    # find an even 2d decomposition that has the most square blocks
    # for a given total number of processors, nproc, such that
    #   nx can be 1 to nproc
    #   nx*ny = nproc
    #   mod(nlats,ny) = mod(nlons,nx) = 0
    #   nscore is the minimum value where
    #     nscore = tmp * tmp where
    #       tmp is (bsize_x/bsize_y - 1)
    #   we want bsize_x/bsize_y to be closest to 1 so subtract
    #   1 and "square it" to create a function that maximizes
    #   squareness when it's minimum.
    # found indicates a decomp has been found, but need to continue
    # to search to find the best decomp.

    if ($found == 0) {
	$nn = 0;
	do {
	    $nn = $nn + 1;
	    $ny = $nn;
	    $nx = int($nprocs/$ny);
	    if ($ny * $nx == $nprocs &&
		$nlats % $ny == 0 &&
		$nlons % $nx == 0) {
		
		$tmp = ($nlons/$nx * $ny/$nlats) - 1.0 ;
		$nscore = $tmp * $tmp;
		if ($nscore < $bscore) {
		    $found = 1;
		    $bscore = $nscore;
		    $mblocks = 1;
		    $dtype = "cartesian";
		    $bsizex = int( $nlons / $nx );
		    $bsizey = int( $nlats / $ny );
		    $nxblocks = $nx;
		    $nyblocks = $ny;
		}
	    }
	    # print "debug $nn $nx $ny $bsizex $bsizey $nscore $bscore $found \n";
	} until ($nn == $nprocs);
    }

    # spacecurve or blockone decomp
    if ($found == 0) 
    {
	my $nscore = 0.0 ;
	my $scok = 0;
	my $fac;
	my $nbx;
	my $nby;
	my $bscore = $nlons * $nlats * $nprocs * 1000;
	my $blksize = ($nlons * $nlats) / ($nprocs);
	$nn = floor(sqrt($blksize)) + 1;
	my $nxnyfact = $nn * 4;
	
	if ($nlats == 1) 
	{
	    $nxnyfact = $blksize;
	}

	$nx = 0;
	do {
	    $nx = $nx + 1;
	    if ($nlons % $nx == 0) {
		$ny = 0;
		do 
		{
		    $ny = $ny + 1;
		    if ($nlats % $ny == 0) {
			$scok = 0;
			# make sure nbx and nby are factorably by only 2, 3, 5 for spacecurve
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
			if ($nbx == 1 && $nby == 1) {
			    $scok = 1;
			}

			# tcraig, somewhat arbitrary scoring system, best is min score
			# min aspect ratio, match blksize, scok good, about 1 block per pe
			my $nblocks = ($nlats * $nlons)/($nx * $ny * $nprocs);
			$nscore = 0.25 * ($nx/$ny + $ny/$nx + ($nx*$ny)/$blksize + $blksize/($nx*$ny))  + 0.5 * (1.0 - $scok) + 0.5 * (1.0/($nblocks) + $nblocks);
			#print "tcxf3 $nlons $nlats $nprocs $blksize $scok $nblocks $nx $ny $nscore \n";
			if ($nscore < $bscore) {
			    $found  = 1;
			    $mblocks = floor(($nlats * $nlons + $nx * $ny * $nprocs - 1)/($nx * $ny * $nprocs));
			    $bscore = $nscore;
			    $bsizex = $nx;
			    $bsizey = $ny;
			    if ($scok == 1) {
				$dtype = "spacecurve";
			    } else {
				$dtype = "blockone";
			    }
			    #print "tcxf4 $bsizex $bsizey $dtype $mblocks \n";
			}
		    }
		} until ($ny > $nxnyfact);
	    }
	}  until ($nx >= $nxnyfact);
    }
    
    if ($found == 1) 
    {
	$decomp{'nlats'}      = $nlats;
	$decomp{'nlons'}      = $nlons;
	$decomp{'maxblocks'}  = $mblocks;
	$decomp{'decomptype'} = $dtype;
	$decomp{'bsize_x'}    = $bsizex;
	$decomp{'bsize_y'}    = $bsizey;
	$decomp{'nx_blocks'}  = $nxblocks;
	$decomp{'ny_blocks'}  = $nyblocks;
    }

    return(%decomp);
}
