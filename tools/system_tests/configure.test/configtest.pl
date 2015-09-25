#!/usr/bin/env perl
use strict;
use warnings;
use XML::LibXML;
use Getopt::Long;
use File::Path;
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
     -generate            generate new baselines in this directory.
     -compare             compare against baselines in this directory.
     -help [or -h]        Print usage to STDOUT (optional).
     -silent [or -s]      Turns on silent mode - only fatal messages issued (optional).
     -verbose [or -v]     Turn on verbose echoing of settings made by create_newcase (optional).
     
EXAMPLES

  ./configtest.pl -cimeroot /path/to/cime -compare /path/to/baselines/cime3.0.1

EOF
}
my ($model, $generate, $compare, $cimeroot, $machdir, $verbose, $silent, $help);
my $logger;
my $loglevel = "INFO";
my $output_dir = getcwd();
sub options{
    GetOptions(
	"cimeroot=s"                => \$cimeroot,
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
    
    if(!defined $generate && !defined $compare){
	die "At least one of the arguments generate and compare must be provided";
    }
    if(defined $compare){
	if(-d $compare){
	    $logger->info( "Compare to baseline directory $compare");
	}else{
	    $logger->logdie ("Could not find baseline directory $compare");
	}
    }
    if(defined $generate){
	if(-d $generate){
	    $logger->warn ("directory $generate already exists");
	}else{
	    mkdir $generate or $logger->logdie ("Could not create directory $generate");
	}
    }else{
	$generate = ".";
    }
    if(defined $output_dir){
	if(-d $output_dir){
	    $logger->warn("Output directory $output_dir already exists");
	}else{
	    mkdir $output_dir or $logger->logdie ("Could not create directory $generate");
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
    my @child = $node->findnodes(".//COMPILERS");
    my @compilers = split(',',$child[0]->textContent());
    
    foreach my $compiler (@compilers){
	    $mach = "yellowstone";
	    $compiler = "intel";
	foreach my $format (qw(make cmake)){
	    
	    my $test_dir = "$output_dir/$mach/$compiler/$format";
	    unless(-d $test_dir){
		mkpath $test_dir or $logger->logdie("Could not create directory $test_dir");
	    }
	    qx( $cimeroot/tools/configure -cimeroot $cimeroot -mach $mach -compiler $compiler -output_dir $test_dir -output_format $format -loglevel $loglevel);
	}
	    last;
    }
    last;
}

