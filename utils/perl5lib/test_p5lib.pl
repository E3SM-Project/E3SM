#!/usr/bin/env perl
# -*- mode: cperl; indent-tabs-mode: nil; cperl-indent-level: 2; -*-
#-----------------------------------------------------------------------------------------------

require 5;

use strict;
use warnings;
use diagnostics;

use Test::More;

BEGIN {
  #
  # setup paths for cime and local utility modules at compile time.
  #
  # Assumes that we are running from cime/perl5lib....
  #
  use Cwd qw(getcwd abs_path);

  my $cwd        = getcwd();
  my $cimeroot   = abs_path("$cwd/../..");

  my @dirs = ("$cwd",
              "$cwd/CPAN",
              "$cimeroot/scripts");

  unshift @INC, @dirs;

  # check for a couple of modules from the paths we added.
  use_ok("Test::Class");
  use_ok("Test::Exception");
}


# ccsm perl modules
require Build::Config;
require Build::NamelistDefinition;
require Build::NamelistDefaults;
require Build::Namelist;
require Streams::TemplateGeneric;
require Tools::SetupTools;

# local perl modules
use Test::Class;
use Test::Exception;

use Test::Class::Load "./t";

Test::Class->runtests;

