#!/usr/bin/env perl 
package test_ModuleLoader;

use Data::Dumper;
use Test::More;
use Test::Exception;
use Module::ModuleLoader;

use parent qw(Test::Class);

sub set_machine
{
	my $self = shift;
	my $machine = shift;
	
	$self->{machine} = $machine;
}
#==============================================================================
# Common test fixtures for all tests. these get reused across all test instances, 
# They MUST be read-only  
#==============================================================================
sub startup: Test(startup => 0) {
	my $self = shift;
	# set up non-machine-specific stuff here. 
	#$self->{scriptsroot} = Cwd::abs_path("../../../scripts/");
	$self->{cimeroot} = Cwd::abs_path("../../../");
	#TODO defined the caseroot?
	
	#if($self->{machine} eq 'goldbach')
	#{
	#	$self->{init_path} = "/usr/share/Modules/init/";
	#	$self->{cmd_path} = "/usr/bin/modulecmd";
	#}
}

#==============================================================================
# 
#==============================================================================
sub shutdown: Test(shutdown) {
	my $self = shift;
}

#==============================================================================
# These get created for each test..
#==============================================================================
sub setup : Test(setup => 0) {
	my $self = shift;
}

#==============================================================================
#
#==============================================================================
sub teardown : Test(setup => 0) {
	my $self = shift;
}


#==============================================================================
# Tests  
#==============================================================================
sub test_new() : Test(1)
{
	my $self = shift;
	#my $machine = $self->{machine};
	my $machine = "goldbach";
	my $compiler = "intel";
	my $mpilib = "openmpi";
	#my $scriptsroot = "../../../scripts";
	#my $cimeroot = "../../../";
	my $cimeroot = "../../";
	
	my $moduleloader = Module::ModuleLoader->new(machine => $machine, compiler => $compiler, mpilib => $mpilib, 
	                                     debug => "FALSE", cimeroot => $cimeroot, caseroot => "$cimeroot/Tools");
	
	isa_ok($moduleloader, "Module::ModuleLoader");
                                         
}

#sub test_moduleInit_brutus() : Test(2)
#{
#	my $self = shift;
#	my $moduleloaderbrutus  = Module::ModuleLoader->new(machine => 'brutus', compiler => 'pgi', mpilib => 'openmpi', 
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#	
#	$moduleloaderbrutus->moduleInit();
#	
#	ok($moduleloaderbrutus->{initpath} eq '/etc/profile.d/modules.perl') || diag($moduleloaderbrutus->{initpath});
#	ok($moduleloaderbrutus->{cmdpath} eq '/usr/bin/modulecmd perl') || diag($moduleloaderbrutus->{cmdpath});
#}
#sub findModulesFromMachinesDir_brutus() : Test(2)
#{
#	my $self = shift;
#	my @expectedbrutusmodules = ( { action => 'purge', actupon => '', seqnum => 1},
#                                  { action => 'load',  actupon => 'intel/10.1.018', seqnum => 2},
#                                  { action => 'load',  actupon => 'mvapich2/1.4rc2', seqnum => 3},
#                                  { action => 'load',  actupon => 'netcdf/4.0.1',    seqnum => 4} );
#
#    my $moduleloaderbrutus  = Module::ModuleLoader->new(machine => 'brutus', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#	
#                                 
#    my @actualmodules = $moduleloaderbrutus->findModulesFromMachinesDir();
#	is_deeply(\@expectedbrutusmodules, \@actualmodules); 
#}
#
#sub test_writeXMLFileForCase_brutus() : Test(1)
#{
#	my $self = shift;
#	
#	my $moduleloader = Module::ModuleLoader->new(machine => 'brutus', compiler => 'intel', mpilib => 'mpich',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#	$moduleloader->writeXMLFileForCase();
#	
#	my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.brutus.xml";
#	open my $EXPECTED, "<", $expectedfile or die "could not open $expectedfile";
#	binmode $EXPECTED;
#	my $expected = do {local $/; <$EXPECTED> };
#	close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#	open my $ACTUAL, "<", $actualfile or die "could not open $actualfile";
#    binmode $ACTUAL;
#	my $actual = do { local $/; <$ACTUAL> };
#	close $ACTUAL;
#	
#	cmp_ok($actual, 'eq', $expected);
#	unlink $actualfile;
#}
#
#sub test_findModulesForCase_brutus() : Test(1)
#{
#	my $self = shift;
#	
#	my $moduleloader = Module::ModuleLoader->new(machine => 'brutus', compiler => 'intel', mpilib => 'mpich',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#	
#	$moduleloader->writeXMLFileForCase();
#	my @actualmodules = $moduleloader->findModulesForCase();
#	
#    my @expectedmodules = ( { action => 'purge', actupon => '', seqnum => 1},
#                                  { action => 'load',  actupon => 'intel/10.1.018', seqnum => 2},
#                                  { action => 'load',  actupon => 'mvapich2/1.4rc2', seqnum => 3},
#                                  { action => 'load',  actupon => 'netcdf/4.0.1',    seqnum => 4} );
#    is_deeply(\@expectedmodules, \@actualmodules);
#	
#}
#
#sub test_writeCshModuleFile_brutus() : Test(1):
#{
#	my $self = shift;
#	my $moduleloader = Module::ModuleLoader->new(machine => 'brutus', compiler => 'intel', mpilib => 'mpich',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#	$moduleloader->findModulesForCase();
#	$moduleloader->writeCshModuleFile();
#	
#	my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.brutus.csh";
#	open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#	binmode $EXPECTED; 
#    my $expected = do { local $/; <$EXPECTED> };
#	close $EXPECTED;
#
#	my $actualfile = "./.env_mach_specific.csh";
#	open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#	binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#	ok($actual eq $expected) || diag "wtf";
#	#unlink $actualfile;
#}
#
#sub test_writeShModuleFile_brutus() : Test(1):
#{
#	my $self = shift;
#	
#	my $moduleloader = Module::ModuleLoader->new(machine => 'brutus', compiler => 'intel', mpilib => 'mpich',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#	$moduleloader->findModulesForCase();
#	$moduleloader->writeShModuleFile();
#
#	my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.brutus.sh";
#	open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#	binmode $EXPECTED; 
#    my $expected = do { local $/; <$EXPECTED> };
#	close $EXPECTED;
#
#	my $actualfile = "./.env_mach_specific.sh";
#	open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#	binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#	ok($actual eq $expected);
#	unlink $actualfile;
#}
#
#sub test_moduleInit_babbage() : Test(2)
#{
#	my $self = shift;
#	my $moduleloader = Module::ModuleLoader->new(machine => 'babbage', compiler => 'intel', mpilib => 'impi', 
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#	
#	$moduleloader->moduleInit();
#	
#	ok($moduleloader->{initpath} eq '/usr/share/Modules/init/perl.pm') || diag($moduleloader->{initpath});
#	ok($moduleloader->{cmdpath} eq '/usr/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_babbage() : Test(1):
#{
#	my $self = shift;
#	my @expectedmodules = ( 
#                            { action => 'unload', actupon => 'intel', seqnum => 1},
#                            { action => 'unload', actupon => 'impi', seqnum => 2},
#                            { action => 'unload', actupon => 'hdf5', seqnum => 3},
#                            { action => 'unload', actupon => 'netcdf', seqnum => 4},
#                            { action => 'load', actupon => 'intel/13.1.2', seqnum => 5},
#                            { action => 'load', actupon => 'impi/4.1.1', seqnum => 6},
#                          );
#	my $moduleloader = Module::ModuleLoader->new(machine => 'babbage', compiler => 'intel', mpilib => 'impi', 
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#	$moduleloader->moduleInit();
#	my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#	is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_babbage() : Test(1):
#{
#	my $self = shift;
#	my @expectedmodules = ( 
#                            { action => 'unload', actupon => 'intel', seqnum => 1},
#                            { action => 'unload', actupon => 'impi', seqnum => 2},
#                            { action => 'unload', actupon => 'hdf5', seqnum => 3},
#                            { action => 'unload', actupon => 'netcdf', seqnum => 4},
#                            { action => 'load', actupon => 'intel/13.1.2', seqnum => 5},
#                            { action => 'load', actupon => 'impi/4.1.1', seqnum => 6},
#                          );
#	my $moduleloader = Module::ModuleLoader->new(machine => 'babbage', compiler => 'intel', mpilib => 'impi', 
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#	my @actualmodules = $moduleloader->findModulesForCase();
#	is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_babbage() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbage', compiler => 'intel', mpilib => 'impi',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.babbage.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#   my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_babbage() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbage', compiler => 'intel', mpilib => 'impi',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.babbage.csh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_babbage() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbage', compiler => 'intel', mpilib => 'impi',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.babbage.sh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_moduleInit_babbageKnc() : Test(2)
#{
#	my $self = shift;
#	my $moduleloader = Module::ModuleLoader->new(machine => 'babbageKnc', compiler => 'intel', mpilib => 'impi', 
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#	$moduleloader->moduleInit();
#	
#	ok($moduleloader->{initpath} eq '/usr/share/Modules/init/perl.pm') || diag($moduleloader->{initpath});
#	ok($moduleloader->{cmdpath} eq '/usr/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_babbageKnc() : Test(1):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                            { action => 'unload', actupon => 'intel', seqnum => 1},
#                            { action => 'unload', actupon => 'impi', seqnum => 2},
#                            { action => 'unload', actupon => 'netcdf', seqnum => 3},
#                            { action => 'load', actupon => 'intel/13.1.2', seqnum => 4},
#                            { action => 'load', actupon => 'impi/4.1.1', seqnum => 5},
#                            { action => 'load', actupon => 'cmake', seqnum => 6},
#                            { action => 'load', actupon => 'netcdf/mic-4.1.3', seqnum => 7},
#                            { action => 'load', actupon => 'pnetcdf/mic-1.5.0', seqnum => 8},
#                          );
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbageKnc', compiler => 'intel13', mpilib => 'impi',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#	#print Dumper \@expectedmodules;
#	#print Dumper \@actualmodules;
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_babbageKnc() : Test(1):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                            { action => 'unload', actupon => 'intel', seqnum => 1},
#                            { action => 'unload', actupon => 'impi', seqnum => 2},
#                            { action => 'unload', actupon => 'netcdf', seqnum => 3},
#                            { action => 'load', actupon => 'intel/13.1.2', seqnum => 4},
#                            { action => 'load', actupon => 'impi/4.1.1', seqnum => 5},
#                            { action => 'load', actupon => 'cmake', seqnum => 6},
#                            { action => 'load', actupon => 'netcdf/mic-4.1.3', seqnum => 7},
#                            { action => 'load', actupon => 'pnetcdf/mic-1.5.0', seqnum => 8},
#                          );
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbageKnc', compiler => 'intel13', mpilib => 'impi',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_babbageKnc() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbageKnc', compiler => 'intel13', mpilib => 'impi',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.babbageKnc.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#   my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_babbageKnc() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbageKnc', compiler => 'intel', mpilib => 'impi',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.babbageKnc.csh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_babbageKnc() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'babbageKnc', compiler => 'intel', mpilib => 'impi',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.babbageKnc.sh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_moduleInit_bluewaters() : Test(2)
#{
#	my $self = shift;
#	my $moduleloader = Module::ModuleLoader->new(machine => 'bluewaters', compiler => 'pgi', mpilib => 'mpich', 
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#	
#	$moduleloader->moduleInit();
#	
#	ok($moduleloader->{initpath} eq '/opt/modules/default/init/perl.pm') || diag($moduleloader->{initpath});
#	ok($moduleloader->{cmdpath} eq '/opt/modules/3.2.10.3/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_bluewaters() : Test(3):
#{
#	my $self = shift;
#	my $moduleloader = Module::ModuleLoader->new(machine => 'bluewaters', compiler => 'pgi', mpilib => 'mpich', 
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#	my @expectedmodules = (
#                             { action => 'rm', actupon => 'PrgEnv-pgi', seqnum => 1}, 
#                             { action => 'rm', actupon => 'PrgEnv-cray', seqnum => 2}, 
#                             { action => 'rm', actupon => 'PrgEnv-gnu', seqnum => 3}, 
#                             { action => 'rm', actupon => 'pgi', seqnum => 4}, 
#                             { action => 'rm', actupon => 'cray', seqnum => 5}, 
#                             { action => 'load', actupon => 'PrgEnv-pgi', seqnum => 6}, 
#                             { action => 'switch', actupon => 'pgi pgi/14.2.0', seqnum => 7}, 
#                             { action => 'load', actupon => 'papi/5.3.2', seqnum => 8}, 
#                             { action => 'switch', actupon => 'cray-mpich cray-mpich/7.0.3', seqnum => 9}, 
#                             { action => 'switch', actupon => 'cray-libsci cray-libsci/12.2.0', seqnum => 10}, 
#                             { action => 'load', actupon => 'torque/5.0.1', seqnum => 11}, 
#                             { action => 'load', actupon => 'cray-netcdf-hdf5parallel/4.3.2', seqnum => 12}, 
#                             { action => 'load', actupon => 'cray-parallel-netcdf/1.5.0', seqnum => 13}, 
#                             { action => 'load', actupon => 'cmake', seqnum => 14}, 
#                             { action => 'rm', actupon => 'darshan', seqnum => 15}, 
#                          );
#     $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    #print Dumper \@expectedmodules;
#    #print Dumper \@actualmodules;
#     is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_bluewaters() : Test(3):
#{
#	my $self = shift;
#	
#	my $moduleloader = Module::ModuleLoader->new(machine => 'bluewaters', compiler => 'pgi', mpilib => 'mpich', 
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    
#    $moduleloader->writeXMLFileForCase();
#	
#	my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.bluewaters.xml";
#	open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#	binmode $EXPECTED;
#	my $expected = do { local $/; <$EXPECTED> };
#	close $EXPECTED;
#
#	my $actualfile = "./env_mach_specific.xml";
#	open(my $ACTUAL, "<", $actualfile) or die " could not open $actualfile";
#	binmode $ACTUAL;
#	my $actual = do { local $/; <$ACTUAL> };
#	close $ACTUAL;
#	cmp_ok($actual, 'eq', $expected);
#	unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_bluewaters() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'bluewaters', compiler => 'pgi', mpilib => 'mpich',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.bluewaters.csh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_bluewaters() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'bluewaters', compiler => 'pgi', mpilib => 'mpich',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.bluewaters.sh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#sub test_findEnvVars_bluewaters() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'bluewaters', compiler => 'pgi', mpilib => 'mpich',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'OMP_STACKSIZE' => '64M',
#          'MPICH_ENV_DISPLAY' => '1',
#          'MPICH_PTL_MATCH_OFF' => '1',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#
#sub test_moduleInit_eastwind() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eastwind', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/etc/profile.d/modules.perl') || diag($moduleloader->{initpath});
#    ok($moduleloader->{cmdpath} eq '/usr/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_eastwind() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eastwind', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    my @expectedmodules = (
#                             { action => 'purge', actupon => '', seqnum => 1},
#                             { action => 'load', actupon => 'pgi/11.3', seqnum => 2},
#                             { action => 'load', actupon => 'mpi/mvapich2/1.5.1p1/pgi11.3', seqnum => 3},
#                             { action => 'load', actupon => 'netcdf/4.1.2/pgi', seqnum => 4},
#                          );
#     $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    #print Dumper \@expectedmodules;
#    #print Dumper \@actualmodules;
#     is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_eastwind() : Test(1):
#{
#    my $self = shift;
#
#   my @expectedmodules = (
#                             { action => 'purge', actupon => '', seqnum => 1},
#                             { action => 'load', actupon => 'pgi/11.3', seqnum => 2},
#                             { action => 'load', actupon => 'mpi/mvapich2/1.5.1p1/pgi11.3', seqnum => 3},
#                             { action => 'load', actupon => 'netcdf/4.1.2/pgi', seqnum => 4},
#                          );
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eastwind', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_eastwind() : Test(3):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eastwind', compiler => 'pgi', mpilib => 'mvapich2',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.eastwind.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die " could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    cmp_ok($actual, 'eq', $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_eastwind() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eastwind', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.eastwind.csh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_eastwind() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eastwind', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.eastwind.sh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_moduleInit_edison() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'edison', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/opt/modules/default/init/perl.pm') || diag($moduleloader->{initpath});
#    ok($moduleloader->{cmdpath} eq '/opt/modules/default/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_edison() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'edison', compiler => 'intel', mpilib => 'mpt',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    my @expectedmodules = (
#                             { action => 'rm', actupon => 'PrgEnv-intel', seqnum => 1},
#                             { action => 'rm', actupon => 'PrgEnv-cray', seqnum => 2},
#                             { action => 'rm', actupon => 'PrgEnv-gnu', seqnum => 3},
#                             { action => 'rm', actupon => 'intel', seqnum => 4},
#                             { action => 'rm', actupon => 'cce', seqnum => 5},
#                             { action => 'rm', actupon => 'cray-parallel-netcdf', seqnum => 6},
#                             { action => 'rm', actupon => 'cray-parallel-hdf5', seqnum => 7},
#                             { action => 'rm', actupon => 'pmi', seqnum => 8},
#                             { action => 'rm', actupon => 'cray-libsci', seqnum => 9},
#                             { action => 'rm', actupon => 'cray-mpich2', seqnum => 10},
#                             { action => 'rm', actupon => 'cray-mpich', seqnum => 11},
#                             { action => 'rm', actupon => 'cray-netcdf', seqnum => 12},
#                             { action => 'rm', actupon => 'cray-hdf5', seqnum => 13},
#                             { action => 'rm', actupon => 'cray-netcdf-hdf5parallel', seqnum => 14},
#                             { action => 'rm', actupon => 'craype-sandybridge', seqnum => 15},
#                             { action => 'rm', actupon => 'craype-ivybridge', seqnum => 16},
#                             { action => 'rm', actupon => 'craype', seqnum => 17},
#                             { action => 'load', actupon => 'PrgEnv-intel', seqnum => 18}, 
#                             { action => 'switch', actupon => 'intel intel/15.0.1.133', seqnum => 19}, 
#                             { action => 'rm', actupon => 'cray-libsci', seqnum => 20}, 
#                             { action => 'use', actupon => '/global/project/projectdirs/ccsm1/modulefiles/edison', seqnum => 21}, 
#                             { action => 'load', actupon => 'esmf/6.2.0-defio-mpi-O', seqnum => 22}, 
#                             { action => 'load', actupon => 'papi/5.3.2', seqnum => 23}, 
#                             { action => 'swap', actupon => 'craype craype/2.1.1', seqnum => 24}, 
#                             { action => 'load', actupon => 'craype-ivybridge', seqnum => 25}, 
#                             { action => 'load', actupon => 'cray-mpich/7.1.1', seqnum => 26}, 
#                             { action => 'load', actupon => 'pmi/5.0.6-1.0000.10439.140.2.ari', seqnum => 27}, 
#                             { action => 'load', actupon => 'cray-netcdf-hdf5parallel/4.3.2', seqnum => 28}, 
#                             { action => 'load', actupon => 'cray-hdf5-parallel/1.8.13', seqnum => 29}, 
#                             { action => 'load', actupon => 'cray-parallel-netcdf/1.5.0', seqnum => 30}, 
#                             { action => 'load', actupon => 'perl/5.20.0', seqnum => 31}, 
#                             { action => 'load', actupon => 'cmake/2.8.11.2', seqnum => 32}, 
#                          );
#     $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    #print Dumper \@expectedmodules;
#    #print Dumper \@actualmodules;
#     is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_edison() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'edison', compiler => 'intel', mpilib => 'mpt',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#    my @expectedmodules = (
#                             { action => 'rm', actupon => 'PrgEnv-intel', seqnum => 1},
#                             { action => 'rm', actupon => 'PrgEnv-cray', seqnum => 2},
#                             { action => 'rm', actupon => 'PrgEnv-gnu', seqnum => 3},
#                             { action => 'rm', actupon => 'intel', seqnum => 4},
#                             { action => 'rm', actupon => 'cce', seqnum => 5},
#                             { action => 'rm', actupon => 'cray-parallel-netcdf', seqnum => 6},
#                             { action => 'rm', actupon => 'cray-parallel-hdf5', seqnum => 7},
#                             { action => 'rm', actupon => 'pmi', seqnum => 8},
#                             { action => 'rm', actupon => 'cray-libsci', seqnum => 9},
#                             { action => 'rm', actupon => 'cray-mpich2', seqnum => 10},
#                             { action => 'rm', actupon => 'cray-mpich', seqnum => 11},
#                             { action => 'rm', actupon => 'cray-netcdf', seqnum => 12},
#                             { action => 'rm', actupon => 'cray-hdf5', seqnum => 13},
#                             { action => 'rm', actupon => 'cray-netcdf-hdf5parallel', seqnum => 14},
#                             { action => 'rm', actupon => 'craype-sandybridge', seqnum => 15},
#                             { action => 'rm', actupon => 'craype-ivybridge', seqnum => 16},
#                             { action => 'rm', actupon => 'craype', seqnum => 17},
#                             { action => 'load', actupon => 'PrgEnv-intel', seqnum => 18},
#                             { action => 'switch', actupon => 'intel intel/15.0.1.133', seqnum => 19},
#                             { action => 'rm', actupon => 'cray-libsci', seqnum => 20},
#                             { action => 'use', actupon => '/global/project/projectdirs/ccsm1/modulefiles/edison', seqnum => 21},
#                             { action => 'load', actupon => 'esmf/6.2.0-defio-mpi-O', seqnum => 22},
#                             { action => 'load', actupon => 'papi/5.3.2', seqnum => 23},
#                             { action => 'swap', actupon => 'craype craype/2.1.1', seqnum => 24},
#                             { action => 'load', actupon => 'craype-ivybridge', seqnum => 25},
#                             { action => 'load', actupon => 'cray-mpich/7.1.1', seqnum => 26},
#                             { action => 'load', actupon => 'pmi/5.0.6-1.0000.10439.140.2.ari', seqnum => 27},
#                             { action => 'load', actupon => 'cray-netcdf-hdf5parallel/4.3.2', seqnum => 28},
#                             { action => 'load', actupon => 'cray-hdf5-parallel/1.8.13', seqnum => 29},
#                             { action => 'load', actupon => 'cray-parallel-netcdf/1.5.0', seqnum => 30},
#                             { action => 'load', actupon => 'perl/5.20.0', seqnum => 31},
#                             { action => 'load', actupon => 'cmake/2.8.11.2', seqnum => 32},
#                          );
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#	print Dumper \@expectedmodules;
#	print Dumper \@actualmodules;
#    is_deeply(\@actualmodules, \@expectedmodules);
#}
#
#sub test_writeXMLFileForCase_edison() : Test(3):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'edison', compiler => 'intel', mpilib => 'mpt',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.edison.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die " could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    cmp_ok($actual, 'eq', $expected);
#    unlink $actualfile;
#}
#sub test_writeCshModuleFile_edison() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'edison', compiler => 'intel', mpilib => 'mpt',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.edison.csh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_edison() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'edison', compiler => 'intel', mpilib => 'mpt',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => ".");
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.edison.sh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_edison() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'edison', compiler => 'intel', mpilib => 'mpich',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'MPICH_ENV_DISPLAY' => '1',
#          'MPICH_VERSION_DISPLAY' => '1',
#          'PERL5LIB' => '/global/project/projectdirs/ccsm1/perl5lib/lib/perl5/5.10.0',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#
#
#sub test_moduleInit_erebus() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'erebus', compiler => 'intel', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/glade/apps/opt/lmod/lmod/init/perl.pm') || diag($moduleloader->{initpath});
#    ok($moduleloader->{cmdpath} eq '/glade/apps/opt/lmod/lmod/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_erebus() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'erebus', compiler => 'intel', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    my @expectedmodules = (
#                             { action => 'purge', actupon => '', seqnum => 1},
#                             { action => 'load', actupon => 'ncarenv/0.0', seqnum => 2},
#                             { action => 'load', actupon => 'ncarbinlibs/0.0', seqnum => 3},
#                             { action => 'load', actupon => 'intel/12.1.4', seqnum => 4},
#                             { action => 'load', actupon => 'ncarcompilers/1.0', seqnum => 5},
#                             { action => 'load', actupon => 'netcdf-mpi/4.2', seqnum => 6},
#                             { action => 'load', actupon => 'pnetcdf/1.3.0', seqnum => 7},
#                          );
#     $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    #print Dumper \@expectedmodules;
#    #print Dumper \@actualmodules;
#     is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_erebus() : Test(1):
#{
#    my $self = shift;
#
#    my @expectedmodules = (
#                             { action => 'purge', actupon => '', seqnum => 1},
#                             { action => 'load', actupon => 'ncarenv/0.0', seqnum => 2},
#                             { action => 'load', actupon => 'ncarbinlibs/0.0', seqnum => 3},
#                             { action => 'load', actupon => 'intel/12.1.4', seqnum => 4},
#                             { action => 'load', actupon => 'ncarcompilers/1.0', seqnum => 5},
#                             { action => 'load', actupon => 'netcdf-mpi/4.2', seqnum => 6},
#                             { action => 'load', actupon => 'pnetcdf/1.3.0', seqnum => 7},
#                          );
#    my $moduleloader = Module::ModuleLoader->new(machine => 'erebus', compiler => 'intel', mpilib => 'ibm',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_erebus() : Test(3):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'erebus', compiler => 'intel', mpilib => 'ibm',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.erebus.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die " could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    cmp_ok($actual, 'eq', $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_erebus() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'erebus', compiler => 'intel', mpilib => 'ibm',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.erebus.csh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_erebus() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'erebus', compiler => 'intel', mpilib => 'ibm',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.erebus.sh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_moduleInit_eos() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eos', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/opt/modules/default/init/perl.pm') || diag($moduleloader->{initpath});
#    ok($moduleloader->{cmdpath} eq '/opt/modules/default/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_eos() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eos', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    my @expectedmodules = (
#                             { action => 'rm', actupon => 'intel', seqnum => 1},
#                             { action => 'rm', actupon => 'cray', seqnum => 2},
#                             { action => 'rm', actupon => 'cray-parallel-netcdf', seqnum => 3},
#                             { action => 'rm', actupon => 'cray-libsci', seqnum => 4},
#                             { action => 'rm', actupon => 'cray-netcdf', seqnum => 5},
#                             { action => 'rm', actupon => 'cray-netcdf-hdf5parallel', seqnum => 6},
#                             { action => 'rm', actupon => 'netcdf', seqnum => 7},
#                             { action => 'load', actupon => 'intel/14.0.2.144', seqnum => 8},
#                             { action => 'load', actupon => 'cray-mpich/7.0.4', seqnum => 9},
#                             { action => 'load', actupon => 'cray-netcdf-hdf5parallel/4.3.2', seqnum => 10},
#                             { action => 'load', actupon => 'cray-parallel-netcdf/1.5.0', seqnum => 11},
#                             { action => 'load', actupon => 'cmake/2.8.11.2', seqnum => 12},
#                          );
#     $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    #print Dumper \@expectedmodules;
#    #print Dumper \@actualmodules;
#     is_deeply(\@actualmodules, \@expectedmodules);
#}
#
#sub test_findModulesForCase_eos() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eos', compiler => 'intel', mpilib => 'mpich',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#    my @expectedmodules = (
#                             { action => 'rm', actupon => 'intel', seqnum => 1},
#                             { action => 'rm', actupon => 'cray', seqnum => 2},
#                             { action => 'rm', actupon => 'cray-parallel-netcdf', seqnum => 3},
#                             { action => 'rm', actupon => 'cray-libsci', seqnum => 4},
#                             { action => 'rm', actupon => 'cray-netcdf', seqnum => 5},
#                             { action => 'rm', actupon => 'cray-netcdf-hdf5parallel', seqnum => 6},
#                             { action => 'rm', actupon => 'netcdf', seqnum => 7},
#                             { action => 'load', actupon => 'intel/14.0.2.144', seqnum => 8},
#                             { action => 'load', actupon => 'cray-mpich/7.0.4', seqnum => 9},
#                             { action => 'load', actupon => 'cray-netcdf-hdf5parallel/4.3.2', seqnum => 10},
#                             { action => 'load', actupon => 'cray-parallel-netcdf/1.5.0', seqnum => 11},
#                             { action => 'load', actupon => 'cmake/2.8.11.2', seqnum => 12},
#                          );
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_eos() : Test(3):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eos', compiler => 'intel', mpilib => 'mpich',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.eos.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die " could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    cmp_ok($actual, 'eq', $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_eos() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eos', compiler => 'intel', mpilib => 'mpich',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.eos.csh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_eos() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eos', compiler => 'intel', mpilib => 'mpich',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.eos.sh";
#    open (my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open (my $ACTUAL, "<", $actualfile) or die "could not open $actualfile, $!";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#sub test_findEnvVars_eos() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'eos', compiler => 'intel', mpilib => 'mpich',
#                                                debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'MPICH_ENV_DISPLAY' => '1',
#          'MPICH_VERSION_DISPLAY' => '1',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#
#
#sub test_moduleInit_hobart() : Test(2)
#{
#    my $self = shift;
#    my $moduleloaderintel  = Module::ModuleLoader->new(machine => 'hobart', compiler => 'intel', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloaderintel->moduleInit();
#
#    ok($moduleloaderintel->{initpath} eq '/usr/share/Modules/init/perl.pm') || diag($moduleloaderintel->{initpath});
#    ok($moduleloaderintel->{cmdpath} eq '/usr/bin/modulecmd perl') || diag($moduleloaderintel->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_hobart() : Test(3):
#{
#    my $self = shift;
#    my @expectedintelmodules = ( {action => 'purge', actupon => '' , seqnum => 1},
#                                 { action => 'load', actupon => 'compiler/intel/15.0.2.164', seqnum => 2}, 
#                                 { action => 'unload', actupon => 'mpi/intel/openmpi-1.8.1-qlc', seqnum => 3},
#                                 { action => 'load', actupon => 'mpi/intel/mvapich2-1.8.1-qlc', seqnum => 4},
#                                 { action => 'load', actupon => 'tool/parallel-netcdf/1.6.1/intel/mvapich2', seqnum => 5}
#                               );
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'intel', mpilib => 'mvapich2',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualintelmodules = $moduleloader->findModulesFromMachinesDir();
#    #print Dumper \@actualintelmodules;
#    is_deeply(\@expectedintelmodules, \@actualintelmodules);
#
#    my @expectedpgimodules = ( { action => 'purge', actupon => '' , seqnum => 1 },
#                               { action => 'load', actupon => 'compiler/pgi/15.1', seqnum => 2},
#                               { action => 'unload', actupon => 'mpi/pgi/openmpi-1.8.1-qlc', seqnum => 3},
#                               { action => 'load', actupon => 'mpi/pgi/mvapich2-1.8.1-qlc', seqnum => 4},
#                               { action => 'load', actupon => 'tool/parallel-netcdf/1.6.1/pgi/mvapich2', seqnum => 5}, 
#                               );
#    my $pgimoduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#    $pgimoduleloader->moduleInit();
#    my @actualpgimodules = $pgimoduleloader->findModulesFromMachinesDir();
#    is_deeply(\@expectedpgimodules, \@actualpgimodules);
#	
#	my @expectednagmodules = ( { action => 'purge', actupon => '', seqnum => 1},
#                               { action => 'load', actupon => 'compiler/nag/6.0', seqnum => 2},
#                               { action => 'load', actupon => 'tool/parallel-netcdf/1.6.1/nag/openmpi', seqnum => 3}, 
#                               { action => 'xmlchange', actupon => 'MPILIB=openmpi', seqnum => 4},
#                               );
#	
#	my $nagmoduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'nag', mpilib => 'openmpi', 
#                                                    debug => 'FALSE', cimeroot => "../../", caseroot => '.');
#	$nagmoduleloader->moduleInit();
#	my @actualnagmodules = $nagmoduleloader->findModulesFromMachinesDir();
#	#print Dumper \@expectednagmodules;
#	#print Dumper \@actualnagmodules;
#	
#	is_deeply(\@expectednagmodules, \@actualnagmodules);
#}
#sub test_writeXMLFileForCase_hobart() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'intel', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.hobart.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> } ;
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_findModulesForCase_hobart() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'intel', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#
#    my @expectedmodules = ( {action => 'purge', actupon => '' , seqnum => 1},
#                                 { action => 'load', actupon => 'compiler/intel/15.0.2.164', seqnum => 2}, 
#                                 { action => 'unload', actupon => 'mpi/intel/openmpi-1.8.1-qlc', seqnum => 3},
#                                 { action => 'load', actupon => 'mpi/intel/mvapich2-1.8.1-qlc', seqnum => 4},
#                                 { action => 'load', actupon => 'tool/parallel-netcdf/1.6.1/intel/mvapich2', seqnum => 5}
#                               );
#    is_deeply(\@expectedmodules, \@actualmodules);
#    #is_deeply(\@expectedmodules, @{$moduleloader->{modulestoload}});
#}
#
#sub test_writeCshModuleFile_hobart() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'intel', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.hobart.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_hobart() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'intel', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.hobart.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_hobart() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hobart', compiler => 'intel', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'P4_GLOBMEMSIZE' => '500000000',
#          'NETCDF_DIR' => '$NETCDF_PATH',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#sub test_moduleInit_gaea() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'gaea', compiler => 'pgi', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/opt/modules/default/init/perl.pm') || diag($moduleloaderintel->{initpath});
#    ok($moduleloader->{cmdpath} eq '/opt/modules/default/bin/modulecmd perl') || diag($moduleloaderintel->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_gaea() : Test(3):
#{
#    my $self = shift;
#    my @expectedmodules = ( 
#                                 { action => 'rm', actupon => 'PrgEnv-pgi' , seqnum => 1},
#                                 { action => 'rm', actupon => 'PrgEnv-cray' , seqnum => 2},
#                                 { action => 'rm', actupon => 'PrgEnv-gnu' , seqnum => 3},
#                                 { action => 'rm', actupon => 'pgi' , seqnum => 4},
#                                 { action => 'rm', actupon => 'cray' , seqnum => 5},
#                                 { action => 'load', actupon => 'PrgEnv-pgi' , seqnum => 6},
#                                 { action => 'switch', actupon => 'pgi pgi/12.5.0' , seqnum => 7},
#                                 { action => 'load', actupon => 'torque/4.1.3' , seqnum => 8},
#                                 { action => 'load', actupon => 'netcdf-hdf5parallel/4.2.0' , seqnum => 9},
#                                 { action => 'load', actupon => 'parallel-netcdf/1.2.0' , seqnum => 10},
#                               );
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'gaea', compiler => 'pgi', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_gaea() : Test(1):
#{
#    my $self = shift;
#
#    my @expectedmodules = (
#                                 { action => 'rm', actupon => 'PrgEnv-pgi' , seqnum => 1},
#                                 { action => 'rm', actupon => 'PrgEnv-cray' , seqnum => 2},
#                                 { action => 'rm', actupon => 'PrgEnv-gnu' , seqnum => 3},
#                                 { action => 'rm', actupon => 'pgi' , seqnum => 4},
#                                 { action => 'rm', actupon => 'cray' , seqnum => 5},
#                                 { action => 'load', actupon => 'PrgEnv-pgi' , seqnum => 6},
#                                 { action => 'switch', actupon => 'pgi pgi/12.5.0' , seqnum => 7},
#                                 { action => 'load', actupon => 'torque/4.1.3' , seqnum => 8},
#                                 { action => 'load', actupon => 'netcdf-hdf5parallel/4.2.0' , seqnum => 9},
#                                 { action => 'load', actupon => 'parallel-netcdf/1.2.0' , seqnum => 10},
#                               );
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'gaea', compiler => 'pgi', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_gaea() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'gaea', compiler => 'pgi', mpilib => 'mpich',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.gaea.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> } ;
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_gaea() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'gaea', compiler => 'pgi', mpilib => 'mpich',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.gaea.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_gaea() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'gaea', compiler => 'pgi', mpilib => 'mpich',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.gaea.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_gaea() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'gaea', compiler => 'pgi', mpilib => 'mpich',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'OMP_STACKSIZE' => '64M',
#          'MPICH_ENV_DISPLAY' => '1',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#sub test_moduleInit_hera() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'hera', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/usr/global/tools/dotkit/init.csh') || diag($moduleloaderintel->{initpath});
#    ok($moduleloader->{cmdpath} eq '') || diag($moduleloaderintel->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_hera() : Test(3):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'use -q', actupon => 'pgi-11.1' , seqnum => 1},
#                                 { action => 'use -q', actupon => 'mvapich2-pgi-1.7' , seqnum => 2},
#                                 { action => 'use -q', actupon => 'netcdf-pgi-4.1.3' , seqnum => 3},
#                          );
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'hera', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_hera() : Test(1):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'use -q', actupon => 'pgi-11.1' , seqnum => 1},
#                                 { action => 'use -q', actupon => 'mvapich2-pgi-1.7' , seqnum => 2},
#                                 { action => 'use -q', actupon => 'netcdf-pgi-4.1.3' , seqnum => 3},
#                          );
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'hera', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_hera() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hera', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.hera.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> } ;
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_hera() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hera', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.hera.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_hera() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hera', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.hera.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_hera() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'hera', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'NETCDF' => '/usr/local/tools/netcdf-pgi-4.1.3/',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#sub test_moduleInit_mira() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'mira', compiler => 'ibm', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/etc/profile.d/00softenv.sh') || diag($moduleloader->{initpath});
#    ok($moduleloader->{cmdpath} eq 'soft') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_mira() : Test(3):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'add', actupon => '+mpiwrapper-xl' , seqnum => 1},
#                                 { action => 'add', actupon => '@ibm-compilers-2015-02' , seqnum => 2},
#                                 { action => 'add', actupon => '+cmake' , seqnum => 3},
#                                 { action => 'add', actupon => '+python' , seqnum => 4},
#                          );
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'mira', compiler => 'ibm', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_mira() : Test(1):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'add', actupon => '+mpiwrapper-xl' , seqnum => 1},
#                                 { action => 'add', actupon => '@ibm-compilers-2015-02' , seqnum => 2},
#                                 { action => 'add', actupon => '+cmake' , seqnum => 3},
#                                 { action => 'add', actupon => '+python' , seqnum => 4},
#                          );
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'mira', compiler => 'ibm', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_mira() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'mira', compiler => 'ibm', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.mira.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> } ;
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_mira() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'mira', compiler => 'ibm', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.mira.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_mira() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'mira', compiler => 'ibm', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.mira.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_mira() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'mira', compiler => 'ibm', mpilib => 'ibm',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'MPI_TYPE_MAX' => '$MPI_TYPE_MAX',
#          'OMP_DYNAMIC' => '$OMP_DYNAMIC',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#sub test_moduleInit_olympus() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'olympus', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/share/apps/modules/Modules/3.2.7/init/perl.pm') || diag($moduleloader->{initpath});
#    ok($moduleloader->{cmdpath} eq '/share/apps/modules/Modules/3.2.7/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_olympus() : Test(3):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'purge', actupon => '' , seqnum => 1},
#                                 { action => 'load', actupon => 'precision/i4' , seqnum => 2},
#                                 { action => 'load', actupon => 'pgi/11.8' , seqnum => 3},
#                                 { action => 'load', actupon => 'mvapich2/1.7' , seqnum => 4},
#                                 { action => 'load', actupon => 'netcdf/4.1.3' , seqnum => 5},
#                          );
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'olympus', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_olympus() : Test(1):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'purge', actupon => '' , seqnum => 1},
#                                 { action => 'load', actupon => 'precision/i4' , seqnum => 2},
#                                 { action => 'load', actupon => 'pgi/11.8' , seqnum => 3},
#                                 { action => 'load', actupon => 'mvapich2/1.7' , seqnum => 4},
#                                 { action => 'load', actupon => 'netcdf/4.1.3' , seqnum => 5},
#                          );
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'olympus', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_olympus() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'olympus', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.olympus.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> } ;
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_olympus() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'olympus', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.olympus.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_olympus() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'olympus', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.olympus.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#sub test_findEnvVars_olympus() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'olympus', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#
#sub test_moduleInit_pleiades_has() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'pleiades-has', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/usr/share/modules/init/perl.pm') || diag($moduleloader->{initpath});
#    ok($moduleloader->{cmdpath} eq '/usr/share/modules/bin/modulecmd perl') || diag($moduleloader->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_pleiades_has() : Test(3):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'purge', actupon => '' , seqnum => 1},
#                                 { action => 'load', actupon => 'comp-intel/2015.0.090' , seqnum => 2},
#                                 { action => 'load', actupon => 'mpi-sgi/mpt.2.11r13' , seqnum => 3},
#                                 { action => 'load', actupon => 'netcdf/4.1.3/intel/mpt' , seqnum => 4},
#                                 { action => 'load', actupon => 'cmake/2.8.12.1' , seqnum => 5},
#                                 { action => 'load', actupon => 'nas' , seqnum => 6},
#                          );
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'pleiades-has', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_pleiades_has() : Test(1):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'purge', actupon => '' , seqnum => 1},
#                                 { action => 'load', actupon => 'comp-intel/2015.0.090' , seqnum => 2},
#                                 { action => 'load', actupon => 'mpi-sgi/mpt.2.11r13' , seqnum => 3},
#                                 { action => 'load', actupon => 'netcdf/4.1.3/intel/mpt' , seqnum => 4},
#                                 { action => 'load', actupon => 'cmake/2.8.12.1' , seqnum => 5},
#                                 { action => 'load', actupon => 'nas' , seqnum => 6},
#                          );
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'pleiades-has', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_pleiades_has() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'pleiades-has', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.pleiades-has.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> } ;
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_pleiades_has() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'pleiades-has', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.pleiades-has.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_pleiades_has() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'pleiades-has', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.pleiades-has.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_pleiades_has() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'pleiades-has', compiler => 'intel', mpilib => 'mpich',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'MPI_GROUP_MAX' => '1024',
#          'MPI_TYPE_MAX' => '100000',
#          'KMP_AFFINITY' => 'noverbose,disabled',
#          'KMP_SCHEDULE' => 'static,balanced',
#          'MPI_TYPE_DEPTH' => '10',
#          'PNETCDF_PATH' => '/home1/fvitt/parallel-netcdf-1.3.1',
#          'OMP_DYNAMIC' => 'FALSE',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#
#sub test_moduleInit_sierra() : Test(2)
#{
#    my $self = shift;
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'sierra', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#
#    ok($moduleloader->{initpath} eq '/usr/global/tools/dotkit/init.csh') || diag($moduleloaderintel->{initpath});
#    ok($moduleloader->{cmdpath} eq '') || diag($moduleloaderintel->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_sierra() : Test(3):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'use -q', actupon => 'pgi-11.1' , seqnum => 1},
#                                 { action => 'use -q', actupon => 'mvapich2-pgi-1.7' , seqnum => 2},
#                                 { action => 'use -q', actupon => 'netcdf-pgi-4.1.3' , seqnum => 3},
#                          );
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'sierra', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_findModulesForCase_sierra() : Test(1):
#{
#    my $self = shift;
#    my @expectedmodules = (
#                                 { action => 'use -q', actupon => 'pgi-11.1' , seqnum => 1},
#                                 { action => 'use -q', actupon => 'mvapich2-pgi-1.7' , seqnum => 2},
#                                 { action => 'use -q', actupon => 'netcdf-pgi-4.1.3' , seqnum => 3},
#                          );
#
#    my $moduleloader  = Module::ModuleLoader->new(machine => 'sierra', compiler => 'pgi', mpilib => 'mvapich2',
#                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->writeXMLFileForCase();
#    my @actualmodules = $moduleloader->findModulesForCase();
#    is_deeply(\@expectedmodules, \@actualmodules);
#}
#
#sub test_writeXMLFileForCase_sierra() : Test(3):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'sierra', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.sierra.xml";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    binmode $EXPECTED;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./env_mach_specific.xml";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    binmode $ACTUAL;
#    my $actual = do { local $/; <$ACTUAL> } ;
#    close $actual;
#    cmp_ok($actual,  'eq',  $expected);
#    unlink $actualfile;
#}
#
#sub test_writeCshModuleFile_sierra() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'sierra', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.sierra.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_sierra() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'sierra', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.sierra.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_sierra() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'sierra', compiler => 'pgi', mpilib => 'mvapich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'NETCDF' => '/usr/local/tools/netcdf-pgi-4.1.3/',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
#sub test_moduleInit_yellowstone() : Test(2)
#{
#	my $self = shift;
#	my $moduleloaderys = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2',
#                                                   debug => "FALSE", cimeroot => "../../", caseroot => ".");
#
#	$moduleloaderys->moduleInit();
#	
#	ok($moduleloaderys->{initpath} eq '/glade/apps/opt/lmod/lmod/init/perl') || diag($moduleloaderys->{initpath});
#	ok($moduleloaderys->{cmdpath} eq '/glade/apps/opt/lmod/lmod/libexec/lmod perl') || diag($moduleloaderys->{cmdpath});
#}
#
#sub test_findModulesFromMachinesDir_yellowstone() : Test(3):
#{
#	my $self = shift;
#	my @expectedintelmpichmodules = ( 
#                            { action => 'purge', actupon => '', seqnum => 1},
#                            { action => 'load', actupon => 'ncarenv/1.0', seqnum => 2},
#                            { action => 'load', actupon => 'ncarbinlibs/1.1', seqnum => 3},
#                            { action => 'load', actupon => 'perlmods', seqnum => 4},
#                            { action => 'load', actupon => 'gmake/4.1', seqnum => 5},
#                            { action => 'load', actupon => 'python', seqnum => 6},
#                            { action => 'load', actupon => 'all-python-libs', seqnum => 7},
#                            { action => 'load', actupon => 'intel/15.0.3', seqnum => 8},
#                            { action => 'load', actupon => 'mkl/11.1.2', seqnum => 9},
#                            { action => 'load', actupon => 'trilinos/11.10.2', seqnum => 10},
#                            { action => 'load', actupon => 'esmf', seqnum => 11},
#                            { action => 'load', actupon => 'esmf-6.3.0rp1-defio-mpi-O', seqnum => 12},
#                            { action => 'load', actupon => 'netcdf-mpi/4.3.3.1', seqnum => 13},
#                            { action => 'load', actupon => 'pnetcdf/1.6.0', seqnum => 14},
#                            { action => 'load', actupon => 'ncarcompilers/1.0', seqnum => 15},
#                            { action => 'load', actupon => 'cmake/2.8.10.2', seqnum => 16},
#						  );
#	#print Dumper \@expectedintelmpichmodules;
#	my $moduleloader = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2', 
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#	$moduleloader->moduleInit();
#	my @actualintelmpichmodules = $moduleloader->findModulesFromMachinesDir();
#    #print Dumper \@expectedintelmpichmodules;
#	#print Dumper \@actualintelmpichmodules;
#	#print "expected: ", ref $expectedintelmpichmodules[0], "\n";
#	#print "actual: ", ref $actualintelmpichmodules[0], "\n";
#	is_deeply(\@actualintelmpichmodules, \@expectedintelmpichmodules, "do modules match");
#
#    my @expectedpgimpichmodules = (
#                            { action => 'purge', actupon => '', seqnum => 1},
#                            { action => 'load', actupon => 'ncarenv/1.0', seqnum => 2},
#                            { action => 'load', actupon => 'ncarbinlibs/1.1', seqnum => 3},
#                            { action => 'load', actupon => 'perlmods', seqnum => 4},
#                            { action => 'load', actupon => 'gmake/4.1', seqnum => 5},
#                            { action => 'load', actupon => 'python', seqnum => 6},
#                            { action => 'load', actupon => 'all-python-libs', seqnum => 7},
#                            { action => 'load', actupon => 'pgi/15.1', seqnum => 8},
#                            { action => 'load', actupon => 'netcdf-mpi/4.3.3.1', seqnum => 9},
#                            { action => 'load', actupon => 'pnetcdf/1.6.0', seqnum => 10},
#                            { action => 'load', actupon => 'ncarcompilers/1.0', seqnum => 11},
#                            { action => 'load', actupon => 'cmake/2.8.10.2', seqnum => 12},
#							);
#	my $moduleloaderpgi = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'pgi', mpilib => 'mpich2', 
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#	$moduleloaderpgi->moduleInit();
#	my @actualpgimpichmodules = $moduleloaderpgi->findModulesFromMachinesDir();
#	is_deeply(\@actualpgimpichmodules, \@expectedpgimpichmodules);
#	
#	my @expectedintelmpiserialdebugmodules = (
#                            { action => 'purge', actupon => '', seqnum => 1},
#                            { action => 'load', actupon => 'ncarenv/1.0', seqnum => 2},
#                            { action => 'load', actupon => 'ncarbinlibs/1.1', seqnum => 3},
#                            { action => 'load', actupon => 'perlmods', seqnum => 4},
#                            { action => 'load', actupon => 'gmake/4.1', seqnum => 5},
#                            { action => 'load', actupon => 'python', seqnum => 6},
#                            { action => 'load', actupon => 'all-python-libs', seqnum => 7},
#                            { action => 'load', actupon => 'intel/15.0.3', seqnum => 8},
#                            { action => 'load', actupon => 'mkl/11.1.2', seqnum => 9},
#                            { action => 'load', actupon => 'trilinos/11.10.2', seqnum => 10},
#                            { action => 'load', actupon => 'esmf', seqnum => 11},
#                            { action => 'load', actupon => 'esmf-6.3.0rp1-defio-uni-g', seqnum => 12},
#                            { action => 'load', actupon => 'netcdf/4.3.3.1', seqnum => 13},
#                            { action => 'load', actupon => 'ncarcompilers/1.0', seqnum => 14},
#                            { action => 'load', actupon => 'cmake/2.8.10.2', seqnum => 15},
#	);
#	my $moduleloadermpiserialdebug = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpi-serial', 
#                                                 debug => 'true', cimeroot => "../../", caseroot => '.');
#	$moduleloadermpiserialdebug->moduleInit();
#	my @actualintelmpiserialdebugmodules = $moduleloadermpiserialdebug->findModulesFromMachinesDir();
#	is_deeply(\@actualintelmpiserialdebugmodules, \@expectedintelmpiserialdebugmodules);
#}
#
#sub test_writeXMLFileForCase_yellowstone() : Test(3):
#{
#	my $self = shift;
#	return;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#    $moduleloader->moduleInit();
#	$moduleloader->writeXMLFileForCase();
#	
#	my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.yellowstone.xml";
#	open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#	binmode $EXPECTED;
#	my $expected = do { local $/; <$EXPECTED> };
#	close $EXPECTED;
#
#	my $actualfile = "./env_mach_specific.xml";
#	open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#	binmode $ACTUAL;
#	my $actual = do { local $/; <$ACTUAL> } ;
#	close $actual;
#	cmp_ok($actual,  'eq',  $expected);
#	unlink $actualfile; 
#}
#
#sub test_findModulesForCase_yellowstone() : Test(1):
#{
#	my $self = shift;
#	
#	my $moduleloader = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#	$moduleloader->writeXMLFileForCase();
#	my @actualintelmpichmodules = $moduleloader->findModulesForCase();
#
#    my @expectedintelmpichmodules = (
#                            { action => 'purge', actupon => '', seqnum => 1},
#                            { action => 'load', actupon => 'ncarenv/1.0', seqnum => 2},
#                            { action => 'load', actupon => 'ncarbinlibs/1.1', seqnum => 3},
#                            { action => 'load', actupon => 'perlmods', seqnum => 4},
#                            { action => 'load', actupon => 'gmake/4.1', seqnum => 5},
#                            { action => 'load', actupon => 'python', seqnum => 6},
#                            { action => 'load', actupon => 'all-python-libs', seqnum => 7},
#                            { action => 'load', actupon => 'intel/15.0.3', seqnum => 8},
#                            { action => 'load', actupon => 'mkl/11.1.2', seqnum => 9},
#                            { action => 'load', actupon => 'trilinos/11.10.2', seqnum => 10},
#                            { action => 'load', actupon => 'esmf', seqnum => 11},
#                            { action => 'load', actupon => 'esmf-6.3.0rp1-defio-mpi-O', seqnum => 12},
#                            { action => 'load', actupon => 'netcdf-mpi/4.3.3.1', seqnum => 13},
#                            { action => 'load', actupon => 'pnetcdf/1.6.0', seqnum => 14},
#                            { action => 'load', actupon => 'ncarcompilers/1.0', seqnum => 15},
#                            { action => 'load', actupon => 'cmake/2.8.10.2', seqnum => 16},
#                          );
#
#	#print Dumper \@expectedintelmpichmodules;
#	#print Dumper \@actualintelmpichmodules;
#    is_deeply(\@actualintelmpichmodules, \@expectedintelmpichmodules);
#}
#
#sub test_writeCshModuleFile_yellowstone() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeCshModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.yellowstone.csh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.csh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_writeShModuleFile_yellowstone() : Test(1):
#{
#    my $self = shift;
#
#    my $moduleloader = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    $moduleloader->writeShModuleFile();
#
#    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.yellowstone.sh";
#    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
#    #my $expected = <$EXPECTED>;
#    my $expected = do { local $/; <$EXPECTED> };
#    close $EXPECTED;
#
#    my $actualfile = "./.env_mach_specific.sh";
#    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
#    my $actual = do { local $/; <$ACTUAL> };
#    close $ACTUAL;
#    ok($actual eq $expected);
#    unlink $actualfile;
#}
#
#sub test_findEnvVars_yellowstone() : Test(1):
#{
#    my $self = shift;
#    my $moduleloader = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2',
#                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
#
#    $moduleloader->moduleInit();
#    $moduleloader->writeXMLFileForCase();
#    $moduleloader->findModulesForCase();
#    my %actualenvs = $moduleloader->findEnvVars();
#    my $expectedenvs = {
#          'OMP_STACKSIZE' => '256M',
#          'MP_LABELIO' => 'yes',
#          'MP_INFOLEVEL' => '2',
#          'MP_SHARED_MEMORY' => 'yes',
#          'MP_EUILIB' => 'us',
#          'MP_MPILIB' => '$MPILIB',
#          'MP_STDOUTMODE' => 'unordered',
#          'MP_RC_USE_LMC' => 'yes',
#        };
#    is_deeply(\%actualenvs, $expectedenvs);
#}
####sub test_loadModules_yellowstone()  : Test(1):
####{
####    my $self = shift;
####
####    my $moduleloader = Module::ModuleLoader->new(machine => 'yellowstone', compiler => 'intel', mpilib => 'mpich2',
####                                                 debug => 'false', cimeroot => "../../", caseroot => '.');
####
####    $moduleloader->moduleInit();
####    $moduleloader->writeXMLFileForCase();
####    $moduleloader->findModulesForCase();
####    $moduleloader->loadModules();
####}
###
###
##
sub test_moduleInit_titan() : Test(2)
{
    my $self = shift;
    my $moduleloaderys = Module::ModuleLoader->new(machine => 'titan', compiler => 'pgi', mpilib => 'mpich2',
                                                   debug => "FALSE", cimeroot => "../../", caseroot => ".");

    $moduleloaderys->moduleInit();

    ok($moduleloaderys->{initpath} eq '/opt/modules/default/init/perl') || diag($moduleloaderys->{initpath});
    ok($moduleloaderys->{cmdpath} eq '/opt/modules/default/bin/modulecmd perl') || diag($moduleloaderys->{cmdpath});
}

sub test_findModulesFromMachinesDir_titan() : Test(1):
{
    my $self = shift;
    my @expectedmodules = (
                                 { action => 'rm', actupon => 'PrgEnv-intel' , seqnum => 1},
                                 { action => 'rm', actupon => 'PrgEnv-pgi' , seqnum => 2},
                                 { action => 'rm', actupon => 'PrgEnv-cray' , seqnum => 3},
                                 { action => 'rm', actupon => 'PrgEnv-gnu' , seqnum => 4},
                                 { action => 'rm', actupon => 'PrgEnv-pathscale' , seqnum => 5},
                                 { action => 'rm', actupon => 'intel' , seqnum => 6},
                                 { action => 'rm', actupon => 'pgi' , seqnum => 7},
                                 { action => 'rm', actupon => 'cray' , seqnum => 8},
                                 { action => 'rm', actupon => 'pathscale' , seqnum => 9},
                                 { action => 'rm', actupon => 'parallel-netcdf' , seqnum => 10},
                                 { action => 'rm', actupon => 'netcdf' , seqnum => 11},
                                 { action => 'rm', actupon => 'cmake' , seqnum => 12},
                                 { action => 'rm', actupon => 'cray-mpich' , seqnum => 13},
                                 { action => 'rm', actupon => 'cray-mpich2' , seqnum => 14},
                                 { action => 'rm', actupon => 'cray-libsci' , seqnum => 15},
                                 { action => 'rm', actupon => 'xt-libsci' , seqnum => 16},
                                 { action => 'rm', actupon => 'cray-netcdf' , seqnum => 17},
                                 { action => 'rm', actupon => 'cray-netcdf-hdf5parallel' , seqnum => 18},
                                 { action => 'rm', actupon => 'cray-parallel-netcdf' , seqnum => 19},
                                 { action => 'load', actupon => 'PrgEnv-pgi' , seqnum => 20},
                                 { action => 'switch', actupon => 'pgi pgi/14.2.0' , seqnum => 21},
                                 { action => 'load', actupon => 'cray-mpich/7.0.4' , seqnum => 22},
                                 { action => 'load', actupon => 'cray-libsci/13.0.1' , seqnum => 23},
                                 { action => 'load', actupon => 'esmf/5.2.0rp2' , seqnum => 24},
                                 { action => 'load', actupon => 'cray-netcdf-hdf5parallel/4.3.2' , seqnum => 25},
                                 { action => 'load', actupon => 'cray-parallel-netcdf/1.5.0' , seqnum => 26},
                                 { action => 'load', actupon => 'subversion' , seqnum => 27},
                                 { action => 'load', actupon => 'cmake/2.8.11.2' , seqnum => 28},
                          );
    my $moduleloader  = Module::ModuleLoader->new(machine => 'titan', compiler => 'pgi', mpilib => 'mpich',
                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
    $moduleloader->moduleInit();
    my @actualmodules = $moduleloader->findModulesFromMachinesDir();
    #print "titan actualmodules\n";
    #print Dumper \@actualmodules;
    is_deeply(\@actualmodules, \@expectedmodules);
}

sub test_findModulesForCase_titan() : Test(1):
{
    my $self = shift;

    my $moduleloader  = Module::ModuleLoader->new(machine => 'titan', compiler => 'pgi', mpilib => 'mpich',
                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
    my @expectedmodules = (
                                 { action => 'rm', actupon => 'PrgEnv-intel' , seqnum => 1},
                                 { action => 'rm', actupon => 'PrgEnv-pgi' , seqnum => 2},
                                 { action => 'rm', actupon => 'PrgEnv-cray' , seqnum => 3},
                                 { action => 'rm', actupon => 'PrgEnv-gnu' , seqnum => 4},
                                 { action => 'rm', actupon => 'PrgEnv-pathscale' , seqnum => 5},
                                 { action => 'rm', actupon => 'intel' , seqnum => 6},
                                 { action => 'rm', actupon => 'pgi' , seqnum => 7},
                                 { action => 'rm', actupon => 'cray' , seqnum => 8},
                                 { action => 'rm', actupon => 'pathscale' , seqnum => 9},
                                 { action => 'rm', actupon => 'parallel-netcdf' , seqnum => 10},
                                 { action => 'rm', actupon => 'netcdf' , seqnum => 11},
                                 { action => 'rm', actupon => 'cmake' , seqnum => 12},
                                 { action => 'rm', actupon => 'cray-mpich' , seqnum => 13},
                                 { action => 'rm', actupon => 'cray-mpich2' , seqnum => 14},
                                 { action => 'rm', actupon => 'cray-libsci' , seqnum => 15},
                                 { action => 'rm', actupon => 'xt-libsci' , seqnum => 16},
                                 { action => 'rm', actupon => 'cray-netcdf' , seqnum => 17},
                                 { action => 'rm', actupon => 'cray-netcdf-hdf5parallel' , seqnum => 18},
                                 { action => 'rm', actupon => 'cray-parallel-netcdf' , seqnum => 19},
                                 { action => 'load', actupon => 'PrgEnv-pgi' , seqnum => 20},
                                 { action => 'switch', actupon => 'pgi pgi/14.2.0' , seqnum => 21},
                                 { action => 'load', actupon => 'cray-mpich/7.0.4' , seqnum => 22},
                                 { action => 'load', actupon => 'cray-libsci/13.0.1' , seqnum => 23},
                                 { action => 'load', actupon => 'esmf/5.2.0rp2' , seqnum => 24},
                                 { action => 'load', actupon => 'cray-netcdf-hdf5parallel/4.3.2' , seqnum => 25},
                                 { action => 'load', actupon => 'cray-parallel-netcdf/1.5.0' , seqnum => 26},
                                 { action => 'load', actupon => 'subversion' , seqnum => 27},
                                 { action => 'load', actupon => 'cmake/2.8.11.2' , seqnum => 28},
                          );

    $moduleloader->moduleInit();
    $moduleloader->writeXMLFileForCase();
    my @actualmodules = $moduleloader->findModulesForCase();


    #print Dumper \@expectedintelmpichmodules;
    #    #print Dumper \@actualintelmpichmodules;
    is_deeply(\@actualmodules, \@expectedmodules);
 }

sub test_writeXMLFileForCase_titan() : Test(1):
{
    my $self = shift;
    my $moduleloader  = Module::ModuleLoader->new(machine => 'titan', compiler => 'pgi', mpilib => 'mpich',
                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
    $moduleloader->moduleInit();
    $moduleloader->writeXMLFileForCase();

    my $expectedfile = "./t/mocks_ModuleLoader/mach_specific.titan.xml";
    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
    binmode $EXPECTED;
    my $expected = do { local $/; <$EXPECTED> };
    close $EXPECTED;

    my $actualfile = "./env_mach_specific.xml";
    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
    binmode $ACTUAL;
    my $actual = do { local $/; <$ACTUAL> } ;
    close $actual;
    cmp_ok($actual,  'eq',  $expected);
    #unlink $actualfile;
}

sub test_writeCshModuleFile_titan() : Test(1):
{
    my $self = shift;


    my $moduleloader  = Module::ModuleLoader->new(machine => 'titan', compiler => 'pgi', mpilib => 'mpich',
                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
    $moduleloader->moduleInit();
    $moduleloader->writeXMLFileForCase();
    $moduleloader->findModulesForCase();
    $moduleloader->writeCshModuleFile();

    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.titan.csh";
    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
    my $expected = do { local $/; <$EXPECTED> };
    close $EXPECTED;

    my $actualfile = "./.env_mach_specific.csh";
    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
    my $actual = do { local $/; <$ACTUAL> };
    close $ACTUAL;
    ok($actual eq $expected);
    #unlink $actualfile;
}

sub test_writeShModuleFile_titan() : Test(1):
{
    my $self = shift;


    my $moduleloader  = Module::ModuleLoader->new(machine => 'titan', compiler => 'pgi', mpilib => 'mpich',
                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
    $moduleloader->moduleInit();
    $moduleloader->writeXMLFileForCase();
    $moduleloader->findModulesForCase();
    $moduleloader->writeShModuleFile();

    my $expectedfile = "./t/mocks_ModuleLoader/env_mach_specific.titan.sh";
    open(my $EXPECTED, "<", $expectedfile) or die "could not open $expectedfile";
    my $expected = do { local $/; <$EXPECTED> };
    close $EXPECTED;

    my $actualfile = "./.env_mach_specific.sh";
    open(my $ACTUAL, "<", $actualfile) or die "could not open $actualfile";
    my $actual = do { local $/; <$ACTUAL> };
    close $ACTUAL;
    ok($actual eq $expected);
    #unlink $actualfile;
}
sub test_findEnvVars_titan() : Test(1):
{
    my $self = shift;
    my $moduleloader = Module::ModuleLoader->new(machine => 'titan', compiler => 'pgi', mpilib => 'mpich',
                                                       debug => "FALSE", cimeroot => "../../", caseroot => '.');
    
    $moduleloader->moduleInit();
    $moduleloader->writeXMLFileForCase();
    $moduleloader->findModulesForCase();
    my %actualenvs = $moduleloader->findEnvVars();
    my $expectedenvs = {
          'MPSTKZ' => '64M',
          'MPICH_CPUMASK_DPSPLAY' => '1',
          'MPICH_ENV_DISPLAY' => '1',
          'MPICH_RANK_REORDER_DISPLAY' => '1',
          'OMP_STACKSIZE' => '64M',
          'MPICH_VERSION_DISPLAY' => '1',
          'CRAY_CPU_TARGET' => 'istanbul',
        };
    is_deeply(\%actualenvs, $expectedenvs);
}
1;
	
