package test_XXX;

# Unit tests for function: XXX

use Data::Dumper;
use Test::More;
use Test::Exception;

use parent qw(Test::Class);

#-------------------------------------------------------------------------------
#
# Common test fixture for all tests:
#
#-------------------------------------------------------------------------------
sub startup : Test(startup => 3) {
  my $self = shift;
  # provide common fixture for all tests, only created once at the
  # start of the tests.
  $self->{cfg} = Build::Config->new("t/input/config_cache_clm4_5_test.xml");
  isnt($self->{cfg}, undef, (caller(0))[3] . " : config object created.");

  $self->{definition} = Build::NamelistDefinition->new("t/input/namelist_definition_clm4_5_test.xml");
  isnt($self->{definition}, undef, (caller(0))[3] . " : namelist_definition object created.");

  $self->{defaults} = Build::NamelistDefaults->new("t/input/namelist_defaults_clm4_5_test.xml");
  isnt($self->{defaults}, undef,  (caller(0))[3] . " : namelist_defaults object created.");
}

sub shutdown : Test(shutdown) {
  # cleanup the single instance test fixtures
}

sub setup : Test(setup => 1) {
  my $self = shift;
  # provide common fixture for all tests, create fresh for each test

  $self->{nl} = Build::Namelist->new();
  isnt($self->{nl}, undef, (caller(0))[3] . " : empty namelist object created.");
}

sub teardown : Test(teardown) {
  # clean up after test
}

#-------------------------------------------------------------------------------
#
# tests
#
#-------------------------------------------------------------------------------
sub test_XXX__YYY : Tests {
  my $self = shift;

  my $msg = "Test that the XXX is set correctly for condition YYY.\n";

  use CLMBuildNamelist qw(message);

  my $opts = { XXX => 1, };
  my $nl_flags = { phys => "clm4_5",
                 };

  CLMBuildNamelist::message($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl});
  my $group = $self->{definition}->get_group_name("XXX");
  my $result = $self->{nl}->get_variable_value($group, "XXX");
  isnt($result, 12345) || diag($msg);
}

#-------------------------------------------------------------------------------

sub test_XXX__WWW : Tests {
  my $self = shift;

  my $msg = "Test that the XXX is set correctly for condition WWW.\n";

  use CLMBuildNamelist qw(message);

  my $opts = { XXX => 1, };
  my $nl_flags = { phys => "clm4_5",
                 };

  CLMBuildNamelist::message($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl});
  my $group = $self->{definition}->get_group_name("XXX");
  my $result = $self->{nl}->get_variable_value($group, "XXX");
  is($result, undef) || diag($msg);
}

#-------------------------------------------------------------------------------
sub test_XXX__ZZZ : Tests {
  my $self = shift;

  my $msg = "Test that using XXX under condition ZZZ results in a fatal error.\n";

  use CLMBuildNamelist qw(XXX);

  my $opts = { XXX => 1, };
  my $nl_flags = { phys => "clm4_5",
                   ZZZ => "ZZZ",
                 };

  dies_ok(sub { CLMBuildNamelist::XXX($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}) }) || diag($msg);
}


1;
