package test_setup_cmdl_run_type;

# Unit tests for function: setup_cmdl_run_type

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
sub test_setup_cmdl_run_type__unset : Tests {
  my $self = shift;

  my $msg = "Test that not setting clm_start_type on the command line results in an error.\n";

  use CLMBuildNamelist qw(setup_cmdl_run_type);

  my $opts = { test => 0, };
  my $nl_flags = { phys => "clm4_5",
                 inputdata_rootdir => 0 };

  dies_ok(sub { CLMBuildNamelist::XXX($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}) }) || diag($msg);

}

#-------------------------------------------------------------------------------

sub test_setup_cmdl_run_type__default : Tests {
  my $self = shift;

  my $msg = "Test that setting clm_start_type to 'default' on the command line results in an error.\n";

  use CLMBuildNamelist qw(setup_cmdl_run_type);

  my $opts = { test => 0,
             clm_start_type => "default" };
  my $nl_flags = { phys => "clm4_5",
                 inputdata_rootdir => 0 };


  dies_ok(sub { CLMBuildNamelist::setup_cmdl_run_type($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl}) }) || diag($msg);
}

#-------------------------------------------------------------------------------

sub test_setup_cmdl_run_type__arbitrary_string : Tests {
  my $self = shift;

  my $msg = "Test that the commandline clm_start_type string is set to ".
    "the namelist and nl_flags for any value except 'default'.\n";

  use CLMBuildNamelist qw(setup_cmdl_run_type);

  my $opts = { test => 0,
             clm_start_type => "foo" };
  my $nl_flags = { phys => "clm4_5",
                 inputdata_rootdir => 0 };

  CLMBuildNamelist::setup_cmdl_run_type($opts, $nl_flags, $self->{definition}, $self->{defaults}, $self->{nl});
  my $group = $self->{definition}->get_group_name("clm_start_type");
  my $result = $self->{nl}->get_variable_value($group, "clm_start_type");
  is($result, "'foo'") || diag($msg);
  is($nl_flags->{'clm_start_type'}, "'foo'") || diag($msg);
}

1;
