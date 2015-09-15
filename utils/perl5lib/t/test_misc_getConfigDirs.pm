#!/usr/bin/env perl
# -*- mode: cperl; indent-tabs-mode: nil; cperl-indent-level: 2; -*-
#-----------------------------------------------------------------------------------------------

package test_misc_getConfigDirs;

# Unit tests for function: Misc::getConfigDirs

use Data::Dumper;
use Test::More;
use Test::Exception;

use parent qw(Test::Class);

use Misc::MiscUtils qw(getConfigDirs);

#-------------------------------------------------------------------------------
#
# Common test fixture for all tests:
#
# NOTE(bja, 2015-06) if you put any test assertions into the common
# fixtures, you need to tell the test system about them by updating
# 'startup => 0' or 'setup => 0' to reflect the number of tests.
#
#-------------------------------------------------------------------------------
sub startup : Test(startup => 0) {
  my $self = shift;
  # provide common fixture for all tests, only created once at the
  # start of the tests. Since these objects resused, they MUST
  # be READ ONLY!
}

sub shutdown : Test(shutdown) {
  # cleanup the single instance test fixtures
}

sub setup : Test(setup => 0) {
  my $self = shift;
  # provide common fixture for all tests, create fresh for each test

}

sub teardown : Test(teardown) {
  # clean up after test
}

#-------------------------------------------------------------------------------
#
# tests
#
#-------------------------------------------------------------------------------
sub test_misc_getConfigDirs__homeSecond : Tests {
  my $self = shift;

  my $msg = "Test that the getConfigDirs returns \$HOME/.cime as the second path.\n";

  my $caseroot = './foo/bar/Tools';
  my $machroot = '/abc/def/machines';
  my $expected = qr/\.cime$/;
  my $result = Misc::MiscUtils::getConfigDirs($caseroot, $machroot);
  like(@{$result}[1], $expected) || diag($msg);
}

#-------------------------------------------------------------------------------

sub test_misc_getConfigDirs__caseToolsFirst : Tests {
  my $self = shift;

  my $msg = "Test that the getConfigDirs returns the case tools directory as the first path.\n";

  my $caseroot = './foo/bar/Tools';
  my $machroot = '/abc/def/machines';
  my $expected = qr/\/Tools$/;
  my $result = Misc::MiscUtils::getConfigDirs($caseroot, $machroot);
  like(@{$result}[0], $expected) || diag($msg);
}



1;
