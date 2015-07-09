#!/usr/bin/env perl 
#

require 5;

use strict;
use warnings;
use diagnostics;

use Test::More;
use Getopt::Long;
BEGIN {
	
	use Cwd qw(getcwd abs_path);
    
	my $cwd = getcwd();
	my $cimeroot = abs_path("$cwd/../../");
	
	my @dirs = ("$cwd",
				"$cwd/CPAN",
                "$cimeroot/scripts");

    unshift @INC, @dirs;
	
	use_ok("Test::Class"); 
	use_ok("Test::Exception"); 
}
 
my $machine;
GetOptions("machine=s" => \$machine);

require Module::ModuleLoader;
require Module::ModuleLoader;

use Test::Class;
use Test::Exception;

use Test::Class::Load "./t/test_ModuleLoader.pm";

Test::Class->runtests;
