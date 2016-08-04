#! /usr/bin/env perl

#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------

{{ batchdirectives }}

use POSIX qw(strftime);
use File::Path;
use File::Copy;
use File::Spec;
use File::Basename;
use XML::LibXML;
my $toolsdir = "./Tools";
push(@INC, $toolsdir);
my $perl5lib = "{{ cimeroot }}/cime/utils/perl5lib";
push(@INC, $perl5lib);
require ConfigCase;
require Run::RunChecks;
require ModuleLoader;


#------------------------------------------------------------------------------
# PE Layout Documentation:
#------------------------------------------------------------------------------
{{ pedocumentation }}
# -------------------------------------------------------------------------
# global data needed by the script, stuff like the max number of threads,
# -------------------------------------------------------------------------

sub main
{
	$ENV{'maxthrds'} = 1;

    # First, get the configuration from every xml file.
    my $buildenv = ConfigCase->new("./Tools/config_definition.xml", "./env_build.xml");
    my %config = ConfigCase->getAllResolved();


    # Change to the case root
    chdir($config{'CASEROOT'});

	qx($config{'BATCHSUBMIT'} ./cesm_tseries_generator.py >> tSeriesStatus 2>&1);
}
main(@ARGV) unless caller;
