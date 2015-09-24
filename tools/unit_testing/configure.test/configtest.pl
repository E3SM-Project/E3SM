#!/usr/bin/env perl
use strict;
use warnings;
use XML::LibXML;
use Getopt::Long;
use File::Path;

# Check for the existence of XML::LibXML in whatever perl distribution happens to be in use.
# If not found, print a warning message then exit.
eval {
    require XML::LibXML;
    XML::LibXML->import();
};
if($@)
{
    die("
  The perl module XML::LibXML is needed for XML parsing in the CIME script system.
  Please contact your local systems administrators or IT staff and have them install it for
  you, or install the module locally.  ");

}
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


sub options{
    GetOptions(
	"cimeroot=s"                => \$cimeroot,
	"model=s"                   => \$model,
	"h|help"                    => \$help,
	"machdir=s"                => \$machdir,
	"s|silent"                  => \$silent,
	"v|verbose"                 => \$verbose,
	"generate=s"              => \$generate,
	"compare=s"               => \$compare,
	)  or usage();

# Give usage message.
    usage() if $help;


# Check for unparsed argumentss
    if (@ARGV) {
	print "ERROR: unrecognized arguments: @ARGV\n";
	usage();
    }
    
    if (! defined $model){
	$model = 'cesm';
    }
    
    if(! defined $cimeroot) {
	$cimeroot = abs_path($ENV{CIMEROOT}) if(defined $ENV{CIMEROOT});
    }
    (-d "$cimeroot")  or  die "Cannot find cimeroot directory \"$cimeroot\" " ;
    
    if(! defined $machdir){
	$machdir  = "$cimeroot/cime_config/${model}/machines";
    }
    
    if(!defined $generate && !defined $compare){
	die "At least one of the arguments generate and compare must be provided";
    }
    if(defined $compare){
	if(-d $compare){
	    print "Compare to baseline directory $compare\n" unless($silent);
	}else{
	    die "Could not find baseline directory $compare";
	}
    }
    if(defined $generate){
	if(-d $generate){
	    warn "directory $generate already exists";
	}else{
	    mkdir $generate or die "Could not create directory $generate";
	}

    }else{
	$generate = ".";
    }
}
options();
#-----------------------------------------------------------------------------------------------
# Make sure we can find required perl modules and configuration files.
# Look for them in the directory that contains the configure script.

# Machines definition file.
my $machine_file = 'config_machines.xml';
(-f "$machdir/$machine_file")  or  
    die("Cannot find machine parameters file $machine_file in directory $machdir ");
    

my $xml = XML::LibXML->new( )->parse_file("$machdir/$machine_file");


foreach my $node ($xml->findnodes(".//machine")){
    my $mach = $node->getAttribute("MACH");
    my @child = $node->findnodes(".//COMPILERS");
    my @compilers = split(',',$child[0]->textContent());
    
    foreach my $compiler (@compilers){
	    $mach = "yellowstone";
	    $compiler = "intel";
	foreach my $format (qw(make cmake)){
	    
	    my $output_dir = "$generate/$mach/$compiler/$format";
	    mkpath $output_dir unless(-d $output_dir);

	    qx( $cimeroot/tools/configure -cimeroot $cimeroot -mach $mach -compiler $compiler -output_dir $output_dir -output_format $format);
	}
	last;
    }
    last;

}

