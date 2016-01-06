#!/usr/bin/env perl
# -*- mode: cperl; cperl-indent-level: 4; indent-tabs-mode: nil; -*-
#==============================================================================
# File: Misc.pm
#
# Misc utility code that didn't fit in any other module.
#
#==============================================================================
use strict;
use warnings;

package Misc::MiscUtils;

# built in system packages
use Exporter;

# user installed packages
use XML::LibXML;

my @EXPORT = qw(getConfigDirs getConfigXMLRoot getBatchSystemType);

#==============================================================================
sub getConfigDirs {
    #
    # Return a list of directories to search for configuration xml files.
    #
    # The desired order is specific --> general, [ caseroot, ~/.cime,
    # machines ], to ensure user overrides are properly picked up.
    #
    my ($caseroot, $machroot) = @_;
    my @configdirs = ($caseroot . '/Tools', "$ENV{'HOME'}\/.cime", $machroot);
    return \@configdirs;
}

#==============================================================================
sub getConfigXMLRoot {
    #
    # Reusable function to extract the xml root from the first file
    # that contains a caller specified xpath query. Example use is to
    # find the xml file that contains a specific machine.
    #
    my $machroot = shift;
    my $caseroot = shift;
    my $filename = shift;
    my $xpathmatch = shift;

    my $configdirs = getConfigDirs($caseroot, $machroot);

    my @mach = undef;
    my $root;
    foreach my $dir (@{$configdirs}) {
        my $configfile = $dir . '/' . $filename;
        if (-f $configfile) {
            my $xml = XML::LibXML->new(no_blanks => 1);
            my $config = $xml->parse_file($configfile);
            $root = $config->getDocumentElement();
            @mach = $root->findnodes($xpathmatch);
            if (@mach) {
                last;
            }
        }
    }
    if(! @mach) {
        die "Could not find '$xpathmatch' in any '$filename' in directories: @{$configdirs}";
    }
    return $root;
}

#==============================================================================
sub getBatchSystemType
{
    #
    # Get the name of the batch system for a specified machine.
    #
    # NOTE(bja, 2015-06) this should probably be in a module called
    # "BatchUtils", but it isn't consistent with the current module
    # with that name.
    #
    my $machine = shift;
    my $machroot = shift;
    my $caseroot = shift;

    my $root = getConfigXMLRoot($machroot, $caseroot, 'config_machines.xml',
                                "/config_machines/machine[\@MACH=\'$machine\']");
    my @batchsystems = $root->findnodes("/config_machines/machine[\@MACH=\'$machine\']/batch_system");
    if(! @batchsystems)
    {
        die "Could not determine batch system type for machine $machine";
    }

    my $batchtype = $batchsystems[0]->getAttribute('type');
    return $batchtype;
}


1;
