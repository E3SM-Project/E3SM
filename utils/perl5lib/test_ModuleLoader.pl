#!/usr/bin/env perl 
#

require 5;

use strict;
use warnings;
use diagnostics;

use Test::More;
Begin {
	
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

require Task::TaskMaker;
require Module::ModuleLoader;

use Test::Class;
use Test::Exception;

use Test::Class::Load "./t/test_TaskMaker.pm";

Test::Class->runtests;
