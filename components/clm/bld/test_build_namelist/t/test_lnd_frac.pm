package test_lnd_frac;

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

  $self->{env_xml} = {};
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
sub test_setup_logic_lnd_frac__fail_if_fatmlndfrc_set : Tests {
  my $self = shift;

  my $msg = "Test that opts->lnd_frac and nl->fatmlndfrc can not be set at the same time.\n";

  use CLMBuildNamelist qw(setup_logic_lnd_frac);

  my $opts = { lnd_frac => 1,
              test => 0,
             };

  # NOTE: don't set inputdata_rootdir so we can tell if the die comes from add_default
  my $nl_flags = { phys => "clm4_5",
                 };

  my $group = $self->{definition}->get_group_name("fatmlndfrc");
  $self->{nl}->set_variable_value($group, "fatmlndfrc", 0);

  dies_ok(sub {CLMBuildNamelist::setup_logic_lnd_frac($opts, $nl_flags, $self->{definition},
                                                      $self->{defaults}, $self->{nl},
                                                      $self->{env_xml});}) || diag($msg);
}

#-------------------------------------------------------------------------------
sub test_setup_logic_lnd_frac__set_fatmlndfrc : Tests {
  my $self = shift;

  my $msg = "Test that fatmlndfrc is set from the stream value supplied by\n" .
    "lnd_frac command line option.\n";

  use CLMBuildNamelist qw(setup_logic_lnd_frac);

  my $opts = { lnd_frac => "dummy_file",
              test => 0,
             };
  my $nl_flags = { phys => "clm4_5",
                   inputdata_rootdir => "/dummy/root/path",
                 };

  CLMBuildNamelist::setup_logic_lnd_frac($opts, $nl_flags, $self->{definition}, $self->{defaults},
                                         $self->{nl}, $self->{env_xml});
  my $group = $self->{definition}->get_group_name("fatmlndfrc");
  my $result = $self->{nl}->get_variable_value($group, "fatmlndfrc");

  is($result, "\'/dummy/root/path/dummy_file\'") || diag($msg);
}

1;
