package test_vichydro;

# Unit tests for function: setup_cmdl_vichydro

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

sub test_setup_cmdl_vichydro__clm4_0 : Tests {
  my $self = shift;

  my $msg = "Test that the setting vichydro is a fatal error in clm4_0.\n";

  use CLMBuildNamelist qw(setup_cmdl_vichydro);

  my $opts = { vichydro => 1, };
  my $nl_flags = { phys => "clm4_0",
                 };

  dies_ok(sub { CLMBuildNamelist::setup_cmdl_vichydro($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}) }) || diag($msg);
}

#-------------------------------------------------------------------------------

sub test_setup_cmdl_vichydro__nl_contradicts_cmdl : Tests {
  my $self = shift;

  my $msg = "Test that the setting vichydro when use_vichydro is false is a fatal error.\n";

  use CLMBuildNamelist qw(setup_cmdl_vichydro);

  my $opts = { vichydro => 1, };
  my $nl_flags = { phys => "clm4_5",
                   vichydro => ".false.",
                 };

  my $group = $self->{definition}->get_group_name("use_vichydro");
  $self->{nl}->set_variable_value($group, "use_vichydro", '.false.' );

  dies_ok(sub { CLMBuildNamelist::setup_cmdl_vichydro($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}) }) || diag($msg);
}

#-------------------------------------------------------------------------------

sub test_setup_cmdl_vichydro__set_use_vichydro : Tests {
  my $self = shift;

  my $msg = "Test that the setting vichydro on the commandline sets use_vichydro ".
    "namelist variable to true.\n";

  use CLMBuildNamelist qw(setup_cmdl_vichydro);

  my $opts = { vichydro => 1, };
  my $nl_flags = { phys => "clm4_5",
                 };

  CLMBuildNamelist::setup_cmdl_vichydro($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl});
  my $group = $self->{definition}->get_group_name("use_vichydro");
  my $result = $self->{nl}->get_variable_value($group, "use_vichydro");
  is($result, '.true.') || diag($msg);
}


1;
