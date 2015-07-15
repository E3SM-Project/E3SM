#!/usr/bin/env perl 
package Module::ModuleLoader;
#------------------------------------------------------------------------------
# Perl module intended for setting up the machine environment (loading modules, soft environment, etc)
# Currently specific to the module system, which as of Oct 2014 is every CESM supported machine save
# mira. 
#------------------------------------------------------------------------------
use strict;
use warnings;
use XML::LibXML;
use Data::Dumper;
use Cwd;

#------------------------------------------------------------------------------
# Constructor.  We need the machine and scriptsroot, and optionally the compiler, 
# mpi library, and whether DEBUG is on.   We also want to initialize the module system so this object
# is ready to load modules.  
#------------------------------------------------------------------------------
sub new 
{
	my ($class, %params) = @_;
	if(! defined $params{'machine'})
	{
		die "ModuleLoader requires a machine argument";
	}
	if(! defined $params{'scriptsroot'})
	{
		die "ModuleLoader requires the scriptsroot";
	}
	if(! defined $params{'caseroot'})
	{
		die "ModuleLoader requires the caseroot";
	}
	my $self = {	
		machine => $params{'machine'} || undef,
		compiler => $params{'compiler'} || undef,
		mpilib =>   $params{'mpilib'} || undef,
		debug  =>   $params{'debug'} || undef,
		scriptsroot => $params{'scriptsroot'} || undef,
		caseroot => $params{'caseroot'} || undef,
	};
	if(! defined $params{'debug'} )
	{
		$self->{debug} = 'FALSE';
	}
	$self->{machineroot} = Cwd::abs_path($self->{scriptsroot}) . "/ccsm_utils/Machines/";
	my $toolsroot = $self->{'caseroot'} . "/Tools";
	push(@INC, $toolsroot);
	$self->{'toolsroot'} = $toolsroot;
	require ConfigCase;
	bless $self, $class;
	#$self->moduleInit();
	return $self;
}

#------------------------------------------------------------------------------
# eval the env_mach_specific file for this machine, which should give us the new environment, 
# then insert the new or modified environment variables into %ENV
#------------------------------------------------------------------------------
sub loadModules()
{
	my $self = shift;
	#my $self = {	
	#	machine => $params{'machine'} || undef,
	#	compiler => $params{'compiler'} || undef,
	#	mpilib =>   $params{'mpilib'} || undef,
	#	debug  =>   $params{'debug'} || undef,
	#	scriptsroot => $params{'scriptsroot'} || undef,
	#	caseroot => $params{'caseroot'} || undef,
	#};
	#eval qx($self->{cmdpath} $cmd);

	if(defined $ENV{CIME_MODULES_LOADED}) {return $self};

	my %oldenv = %ENV;
	my $envfile = $self->{'caseroot'} . "/env_mach_specific";
	my $cshenv = "env ";
    $cshenv .= "COMPILER=". $self->{'compiler'} . " " ;
    $cshenv .= "MPILIB=". $self->{'mpilib'} . " " ;
    $cshenv .= "DEBUG=". $self->{'debug'} . " " ;
    $cshenv .= "CASEROOT=". $self->{'caseroot'} . " " ;
    $cshenv .= "PERL=TRUE ";
	#my $cmd = $cshenv . " " . $envfile;
	#my $cmd = $cshenv . " " . $envfile . " && printenv";
	my $cmd = $cshenv . " " . $envfile ;
	#print "running command $cmd\n";
	my @output;
	eval { @output = qx($cmd);};
	#eval { $out = `$cmd`;};
	#$out = `$cmd`;
	chomp @output;
	#print Dumper \@output;
	if($?)
	{
		die "could not load modules for machine $self->{'machine'}";
	}
	my %newenv;
	
	#foreach my $line(@output)
	for my $i(0 .. $#output)
	{
		chomp $output[$i];
		if($output[$i] !~ /=/)
		{
			$output[$i] = '';
		}
		if($output[$i] =~ /BASH_FUNC/i)
		{
			$output[$i] = '';
		}
	}
	#print Dumper \@output;
	foreach my $line(@output)
	{
		if(length($line) > 0)
		{
	            chomp $line;
	            #print "line: $line\n";
	            my ($key, $value) = split('=', $line, 2);
	            $newenv{$key} = $value;
                }
	}
	my %newbuildenv;
	
	foreach my $k(keys %newenv)
	{
		if(! defined $oldenv{$k})
		{
		    #print "new env var: $newenv{$k}\n";
			$newbuildenv{$k} = $newenv{$k};
			$ENV{$k} = $newenv{$k};
		}
		if(defined $oldenv{$k} && $newenv{$k} ne $oldenv{$k})
		{
		    #print "modified env var: $newenv{$k}\n";
			$newbuildenv{$k} = $newenv{$k};
			$ENV{$k} = $newenv{$k};
		}
	}
	$ENV{CIME_MODULES_LOADED} = 1;
	return %newbuildenv;

}
1;
