#!/usr/bin/env perl

# Test methods of the Config objects.

#########################

use Test::More tests => 27;  # use "no_plan" until number of tests is determined

#########################

use strict;

use lib "..";
BEGIN {use_ok('Build::Config')};

# check that XML::Lite is being found
BEGIN {use_ok('XML::Lite')};

# create Build::Config object
my $cfg = Build::Config->new("config_definition.xml", "config_setup_eul.xml");
isa_ok($cfg, "Build::Config", "created configuration object");

# check entries in definition file
is($cfg->get('cam_bld'), ".", 'checking config_definition default values');
is($cfg->get('usr_src'), "", 'checking config_definition default values');
is($cfg->get('cpl'), "0", 'checking config_definition default values');
is($cfg->get('phys'), "cam", 'checking config_definition default values');
is($cfg->get('chem'), "none", 'checking config_definition default values');
is($cfg->get('ocn'), "dom", 'checking config_definition default values');
is($cfg->get('target_os'), "", 'checking config_definition default values');

# check trapping of parameters not in the definition file
eval { $cfg->get('xxx'); };
like( $@, qr/ERROR: unknown parameter name: xxx/, 'get() trapped unkown parameter');

# check entries in setup file
is($cfg->get('dyn'), "eul", 'checking config_setup values');
is($cfg->get('res'), "64x128", 'checking config_setup values');

# test setting values
$cfg->set('dyn', 'fv');
$cfg->set('res', '1.9x2.5');
is($cfg->get('dyn'), "fv", 'check setting values');
is($cfg->get('res'), "1.9x2.5", 'check setting values');

# test setting numeric values
my $val = 0;
$cfg->set('cpl', $val);
is($cfg->get('cpl'), "0", 'check setting numeric values');

# Check getting list of valid values
my @expect = ("none","waccm_mozart","waccm_ghg","trop_mozart");
is( $cfg->get_valid_values("chem" ), @expect, "check that chem valid values match" );
my @expect;
is( $cfg->get_valid_values("cam_bld" ), @expect, "check that cam_bld valid values match" );
my @expect = (0,1);
is( $cfg->get_valid_values("cpl" ), @expect, "check that cpl valid values match" );
my @expect = ("one_value");
is( $cfg->get_valid_values("test_one_valid_value" ), @expect, 
    "check that test_one_value valid values match" );

# test setting list of values
$cfg->set('prog_aero', 'dust,seasalt,sulfur');
is($cfg->get('prog_aero'), "dust,seasalt,sulfur", 'check setting list values');

# test setting illegal value in list
eval { $cfg->set('prog_aero', 'caer,caer4'); };
like( $@, qr/ERROR: caer,caer4 is not a valid value for parameter prog_aero:/,
      'set() trapped illegal list value');

# check that set() traps invalid parameter names and values
eval { $cfg->set('xxx', ""); };
like( $@, qr/ERROR: parameter name xxx is not in the configuration definition/,
      'set() trapped unkown parameter name');

eval { $cfg->set('dyn', "xxx"); };
like( $@, qr/ERROR: xxx is not a valid value for parameter dyn: valid values are eul,sld,fv,homme/,
      'set() trapped invalid parameter value');

eval { $cfg->set('dyn', "eul,sld"); };
like( $@, qr/ERROR: eul,sld is not a valid value for parameter dyn: valid values are eul,sld,fv,homme/,
      'set() trapped invalid list parameter value');

# check write_file method -- use new Build::Config object so can compare input
# definition file with output file
# *** Note that the config_definition.xml file used for testing is specially formatted to
#     match the output of the write_file method.  The entry elements are alphabetized, 
#     there are valid_values attributes even when they contain no values, and all excess
#     whitespace is eliminated from the start tags.  These restrictions are not necessary
#     in configure definition xml files.
$cfg = Build::Config->new("config_definition.xml");
$cfg->write_file("config_cache.xml");
my $rval = `diff config_definition.xml config_cache.xml`;
is($rval, "", 'check writing output file');

# create a config object without reading in files and make sure it fails

eval { my $cfg = Build::Config->new(); };
like( $@, qr/ERROR: Build::Config new method requires a definition file/,
      'check that can`t new a config object without a definition file');

print "\nSuccessfully ran all tests\n";
