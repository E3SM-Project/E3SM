#!/usr/bin/env perl
# -*- mode: cperl; cperl-indent-level: 4; indent-tabs-mode: nil; -*-
#==============================================================================
# File: UserdefinedUtils.pm
#
# Purpose: Provide a simpler way to configure userdefined machines
# without modifying the cime machines directory (don't want everyones
# machine included in cime) or hand modifying xml files in case
# directories (not possible for automated testing).
#
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
    my ($caseroot, $machroot) = @_;
    my @configdirs = ( "$ENV{'HOME'}\/.cime", $machroot, $caseroot . '/Tools' );
    return \@configdirs;
}

#==============================================================================
sub getConfigXMLRoot {
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
        die "Could not find '$xpathmatch' in any '$filename' in directories: $configdirs";
    }
    return $root;
}

#==============================================================================
sub getBatchSystemType
{
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
