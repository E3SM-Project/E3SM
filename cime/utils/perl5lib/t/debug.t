#!/usr/bin/env perl

# Debug methods of the Config objects.

use strict;

use lib "..";

use Build::Config;

# create Build::Config object
my $cfg = Build::Config->new("config_definition.xml", "config_setup_eul.xml");


# test setting illegal value in list
$cfg->set('prog_aero', 'caer,caer4');
my $list_val = $cfg->get('prog_aero');

print "value= $list_val\n";
