#!/usr/bin/env perl
# -*- mode: cperl; indent-tabs-mode: nil; cperl-indent-level: 2; -*-
#-----------------------------------------------------------------------------------------------

package test_template;

# Unit tests for function: template

use Data::Dumper;
use Test::More;
use Test::Exception;

use parent qw(Test::Class);

#-------------------------------------------------------------------------------
#
# Common test fixture for all tests:
#
# NOTE(bja, 2015-06) if you put any test assertions into the common
# fixtures, you need to tell the test system about them by updating
# 'startup => 0' or 'setup => 0' to reflect the number of tests.
#
#-------------------------------------------------------------------------------
sub startup : Test(startup => 1) {
  my $self = shift;
  # provide common fixture for all tests, only created once at the
  # start of the tests. Since these objects resused, they MUST
  # be READ ONLY!

  my $msg = "Test that some startup code ran correctly.\n";
  my $expected = 1.2345;
  my $result = 3.14159;
  isnt($result, $expected) || diag($msg);
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
sub test_template__YYY : Tests {
  my $self = shift;

  my $msg = "Test that the template is set correctly for condition YYY.\n";

  my $expected = "abc";
  my $result = "def";
  isnt($result, $expected) || diag($msg);
}

#-------------------------------------------------------------------------------

sub test_template__WWW : Tests {
  my $self = shift;

  my $msg = "Test that the template is set correctly for condition WWW.\n";

  my $result = undef;
  my $expected = undef;
  is($result, $expected) || diag($msg);
}

#-------------------------------------------------------------------------------
sub test_template__ZZZ : Tests {
  my $self = shift;

  my $msg = "Test that using template under condition ZZZ results in a fatal error.\n";


  dies_ok(sub { die; }) || diag($msg);
}


1;
