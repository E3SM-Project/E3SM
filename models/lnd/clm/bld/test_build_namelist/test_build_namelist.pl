#!/usr/bin/env perl
# -*- mode: cperl -*-
#-----------------------------------------------------------------------------------------------

require 5;

use strict;
use warnings;
use diagnostics;

use Test::More;

BEGIN {
  #
  # setup paths for cesm and local utility modules at compile time.
  #
  # Assumes that we are running from clm/bld/test....
  #
  use Cwd qw(getcwd abs_path);

  my $cwd = getcwd();
  my $cesm_dir;
  my $clm_dir;

  if ($cwd =~ m|((/.+)/models/lnd/clm)/bld|) {
    $cesm_dir = $2;
    $clm_dir = $1;
  }

  my $cesm_tools = "$cesm_dir/scripts/ccsm_utils/Tools/perl5lib";
  my $cpan_dir = "$cwd/perl5lib";
  my @dirs = ("../", $cpan_dir, $cesm_tools);

  unshift @INC, @dirs;

  # check for a couple of modules from the paths we added.
  use_ok("XML::Lite");
  use_ok("CLMBuildNamelist");
  use_ok("Test::Class");
}


# ccsm perl modules
require XML::Lite;
require Build::Config;
require Build::NamelistDefinition;
require Build::NamelistDefaults;
require Build::Namelist;
require Streams::TemplateGeneric;

# local perl modules
use Test::Class;
use Test::Exception;

require CLMBuildNamelist;

use Test::Class::Load "./t";

Test::Class->runtests;

