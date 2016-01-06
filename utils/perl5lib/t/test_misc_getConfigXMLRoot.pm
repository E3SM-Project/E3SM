#!/usr/bin/env perl
# -*- mode: cperl; indent-tabs-mode: nil; cperl-indent-level: 2; -*-
#-----------------------------------------------------------------------------------------------

package test_misc_getConfigXMLRoot;

# Unit tests for function: XXX

use Data::Dumper;
use Test::More;
use Test::Exception;

use parent qw(Test::Class);

use Misc::MiscUtils qw(getConfigXMLRoot);

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
  use File::Path qw(make_path);
  use Cwd qw(getcwd abs_path);

  my $cwd = getcwd();
  $self->{'tmp_dir'} = "$cwd/tmp";
  $self->{'caseroot'} = "$self->{'tmp_dir'}/caseroot";
  make_path("$self->{'caseroot'}/Tools");
  my $filename = "$self->{'caseroot'}/Tools/tmp_mach.xml";
  my $tmp_config = << "MSG";
<?xml version="1.0"?>
<config_machines>
    <machine MACH="foo">
    </machine>
    <machine MACH="bar">
    </machine>
    <machine MACH="zombie">
        <action>shamble quickly</action>
    </machine>
</config_machines>
MSG
  open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
  print $fh $tmp_config;
  close $fh;

  $self->{'machroot'} = "$self->{'tmp_dir'}/machroot";
  make_path($self->{'machroot'});

  $filename = "$self->{'machroot'}/tmp_mach.xml";
  $tmp_config = << "MSG";
<?xml version="1.0"?>
<config_machines>
    <machine MACH="baz">
    </machine>
    <machine MACH="bob">
    </machine>
    <machine MACH="zombie">
        <action>eat brains</action>
    </machine>
</config_machines>
MSG
  open($fh, '>', $filename) or die "Could not open file '$filename' $!";
  print $fh $tmp_config;
  close $fh;
}

sub shutdown : Test(shutdown) {
  # cleanup the single instance test fixtures
  my $self = shift;

  use File::Path qw(remove_tree);
  remove_tree($self->{'tmp_dir'});
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
sub test_misc_getConfigXMLRoot__returnsValidElementCaseRoot : Tests {
  my $self = shift;

  my $msg = "Test that the xml root element is returned when the element is found only in a caseroot file.\n";

  my $expected = undef;
  my $filename = "tmp_mach.xml";
  my $xpathmatch = "/config_machines/machine[\@MACH=\"foo\"]";
  my $result = Misc::MiscUtils::getConfigXMLRoot($self->{'machroot'},
                                                 $self->{'caseroot'},
                                                 $filename,
                                                 $xpathmatch);
  isnt($result, $expected) || diag($msg);
}

sub test_misc_getConfigXMLRoot__returnsValidElementMachRoot : Tests {
  my $self = shift;

  my $msg = "Test that the xml root element is returned when the element is found only in a machroot file.\n";

  my $expected = undef;
  my $filename = "tmp_mach.xml";
  my $xpathmatch = "/config_machines/machine[\@MACH=\"baz\"]";
  my $result = Misc::MiscUtils::getConfigXMLRoot($self->{'machroot'},
                                                 $self->{'caseroot'},
                                                 $filename,
                                                 $xpathmatch);
  isnt($result, $expected) || diag($msg);
}

sub test_misc_getConfigXMLRoot__returnsValidElementMachRootCaseRoot : Tests {
  my $self = shift;

  my $msg = "Test that the xml root element from caseroot is returned when the element is in machroot and caseroot.\n";

  my $expected = "shamble quickly";
  my $filename = "tmp_mach.xml";
  my $xpathmatch = "/config_machines/machine[\@MACH=\"zombie\"]";
  my $xmlroot = Misc::MiscUtils::getConfigXMLRoot($self->{'machroot'},
                                                  $self->{'caseroot'},
                                                  $filename,
                                                  $xpathmatch);
  my $result = $xmlroot->findnodes("/config_machines/machine[\@MACH=\"zombie\"]/action");
  is($result, $expected) || diag($msg);
}

sub test_misc_getConfigXMLRoot__returnsInvalid : Tests {
  my $self = shift;

  my $msg = "Test that the xml root element is undefined when the request xpath query doesn't match any file\n.";

  my $filename = "tmp_mach.xml";
  my $xpathmatch = "/some/invalid[\@junk=\'string\']";
  dies_ok(sub {Misc::MiscUtils::getConfigXMLRoot($self->{'machroot'},
                                            $self->{'caseroot'},
                                            $filename,
                                            $xpathmatch) }) || diag($msg);
}


1;
