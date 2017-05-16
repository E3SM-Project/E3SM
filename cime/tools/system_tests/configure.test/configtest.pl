#!/usr/bin/env perl
use strict;
use warnings;
use XML::LibXML;
use Getopt::Long;
use File::Path;
use File::Basename qw(basename dirname);
use Cwd;

#-----------------------------------------------------------------------------------------------
if ($#ARGV == -1) {
    usage();
}

#-----------------------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     configtest.pl [options]
OPTIONS
     User supplied values are denoted in angle brackets (<>).  Any value that contains
     white-space must be quoted.  Long option names may be supplied with either single
     or double leading dashes.  A consequence of this is that single letter options may
     NOT be bundled.

     -cimeroot            Specify the toplevel cime directory.
                          default: use CIMEROOT environment variable
     -mach_dir <path>     Specify the locations of the Machines directory (optional).
     -model    <name>     Specify the cime model on which to conduct tests.
                          default: cesm
     -baselineroot       Root directory to read or write baselines to.
     -compare             compare against baselines in this directory (under baselineroot).
     -generate             write new baselines to this directory (under baselineroot).
     -help [or -h]        Print usage to STDOUT (optional).

EXAMPLES

  ./configtest.pl -cimeroot /path/to/cime -compare cime3.0.1 -baselineroot /path/to/baselines/

EOF
}
my ($model, $baselineroot, $generate, $compare, $cimeroot, $machdir, $verbose, $silent, $help);
my $logger;
my $loglevel = "INFO";
my $output_dir = getcwd();
my @time = localtime();
$output_dir .= "/$time[5]$time[4]$time[3]$time[2]$time[1]";
sub options{
    GetOptions(
	"cimeroot=s"                => \$cimeroot,
        "baselineroot=s"          => \$baselineroot,
	"model=s"                   => \$model,
	"h|help"                    => \$help,
	"machdir=s"                => \$machdir,
	"loglevel=s"               => \$loglevel,
	"generate=s"              => \$generate,
	"compare=s"               => \$compare,
        "output_dir=s"           => \$output_dir,
	)  or usage();


    if(! defined $cimeroot) {
	$cimeroot = abs_path($ENV{CIMEROOT}) if(defined $ENV{CIMEROOT});
    }
    (-d "$cimeroot")  or  die "Cannot find cimeroot directory \"$cimeroot\" " ;

    my @dirs = ("$cimeroot/utils/perl5lib" );

    unshift @INC, @dirs;
    require Log::Log4perl;
    require File::Copy::Recursive;
    require File::DirCompare;

    my $level = Log::Log4perl::Level::to_priority($loglevel);
    Log::Log4perl->easy_init({level=>$level,
			  layout=>'%m%n'});

    $logger = Log::Log4perl::get_logger();
# Give usage message.
    usage() if $help;


# Check for unparsed argumentss
    if (@ARGV) {
	$logger->error( "ERROR: unrecognized arguments: @ARGV");
	usage();
    }

    if (! defined $model){
	$model = 'cesm';
    }

    if(! defined $machdir){
	$machdir  = "$cimeroot/cime_config/${model}/machines";
    }
    if(!defined $baselineroot){
	$logger->logdie("baselineroot argument must be provided");
    }
    if(!defined $generate && !defined $compare){
	$logger->logdie( "At least one of the arguments generate and compare must be provided");
    }
    if(defined $compare){
	if(-d "$baselineroot/$compare"){
	    $logger->info( "Compare to baseline directory $baselineroot/$compare");
	}else{
	    $logger->logdie ("Could not find baseline directory $compare");
	}
    }
    if(defined $generate){
	if(-d "$baselineroot/$generate"){
	    $logger->logdie ("Directory $baselineroot/$generate already exists.");
	}
    }
    if(defined $output_dir){
	if(-d $output_dir){
	    $logger->warn("Output directory $output_dir already exists");
	}else{

	}
    }
}
options();
#-----------------------------------------------------------------------------------------------
# Make sure we can find required perl modules and configuration files.
# Look for them in the directory that contains the configure script.

# Machines definition file.
my $machine_file = 'config_machines.xml';
(-f "$machdir/$machine_file")  or
    $logger->logdie("Cannot find machine parameters file $machine_file in directory $machdir ");


my $xml = XML::LibXML->new( )->parse_file("$machdir/$machine_file");


foreach my $node ($xml->findnodes(".//machine")){
    my $mach = $node->getAttribute("MACH");
    next if ($mach eq "userdefined");
    my @child = $node->findnodes(".//COMPILERS");
    my @compilers = split(',',$child[0]->textContent());

    foreach my $compiler (@compilers){
	foreach my $format (qw(make cmake)){

	    my $test_dir = "$output_dir/$mach/$compiler/$format";
	    unless(-d $test_dir){
		mkpath $test_dir or $logger->logdie("Could not create directory $test_dir");
	    }
	    qx( $cimeroot/tools/configure -cimeroot $cimeroot -mach $mach -compiler $compiler -output_dir $test_dir -output_format $format -loglevel $loglevel);
	}
    }
}

if(defined $generate){
     File::Copy::Recursive::dircopy($output_dir, "$baselineroot/$generate") or $logger->logdie("Could not copy $output_dir to $baselineroot/$generate");
}


if(defined $compare){
    my $result = "PASS";
    my $filediff =0;
    my $nobase=0;
    my $notest = 0;


    $logger->info("\n\nComparing $output_dir to $baselineroot/$compare");
    File::DirCompare->compare("$output_dir","$baselineroot/$compare",
		     sub {
			 my($a, $b) = @_;
			 if( !$b){
			     my $base = basename($a);
			     my $dir = dirname($a);
			     $logger->warn("  File $base only exists in $dir");
			     $notest++;
			 }elsif( !$a){
			     my $base = basename($b);
			     my $dir = dirname($b);
			     $logger->warn("  File $base only exists in $dir");
			     $nobase++;
			 }else{
			     $logger->error("  File contents for $a and $b differ.\n");
			     $filediff++;
			 }
		     });

    if($filediff >0){
	print "CONFIGTEST RESULT: FAIL $notest $nobase $filediff\n";
    }else{
	print "CONFIGTEST RESULT: PASS\n";
    }


}
