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
	if(! defined $params{'cimeroot'})
	{
		die "ModuleLoader requires the cimeroot";
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
		cimeroot => $params{'cimeroot'} || undef,
		caseroot => $params{'caseroot'} || undef,
	};
	if(! defined $params{'debug'} )
	{
		$self->{debug} = 'FALSE';
	}
	$self->{machroot} = Cwd::abs_path($self->{cimeroot}) . "/machines/";
	print "machroot: $self->{machroot}\n";
	my $toolsroot = $self->{'caseroot'} . "/Tools";
	push(@INC, $toolsroot);
	$self->{'toolsroot'} = $toolsroot;
	#require ConfigCase;
	bless $self, $class;
	#$self->moduleInit();
	return $self;
}

#------------------------------------------------------------------------------
# eval the env_mach_specific file for this machine, which should give us the new environment, 
# then insert the new or modified environment variables into %ENV
#------------------------------------------------------------------------------
sub loadModulesCshEval()
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
	
	return %newbuildenv;

}

sub moduleInit()
{
	my $self = shift;
	#my $envfile = $self->{'modulepath'} . 
	my $configmachinesfile = "$self->{machroot}/config_machines.xml";
	my $machine = $self->{machine};
	if(! -e $configmachinesfile)
	{
        die "$configmachinesfile not found\n";     
	}
	else
	{
		$self->{configmachinesfile} = $configmachinesfile;
	}
	
	# Set up the XML::LibXML parser..
	my $parser = XML::LibXML->new(no_blanks => 1);
	#my $xml = $parser->parse_file($$configmachinesfile);
	my $xml = $parser->parse_file($self->{'configmachinesfile'});
	$self->{configmachinesroot} = $xml;
	
	# Get the init_path.  Module systems usually have an 'init' script for 
	# various scripting languages, we need to get this path from config_machines
	
	my @initnodes = $xml->findnodes("//machine[\@MACH=\'$machine\']/module_system/init_path");
	foreach my $initnode(@initnodes)
	{
		$self->{initpath} = $initnode->textContent();
	}
	if(! defined $self->{initpath})
	{
		die "the module init path could not be found for the machine $machine\n";
	}

	my @cmdnodes = $xml->findnodes("//machine[\@MACH=\'$machine\']/module_system/cmd_path");
	foreach my $cmdnode(@cmdnodes)
	{
		$self->{cmdpath} = $cmdnode->textContent();
	}
	if(! defined $self->{cmdpath})
	{
		die "the module cmd path could not be found for the machine $machine\n";
	}

	#$self->{initpath} .= "perl.pm";
	#$self->{cmdpath} .=  " perl";
}

sub findModulesFromMachinesDir()
{
	my $self = shift;
	
	# Get the compiler, mpilib, machine, and debug status. 
	my $compiler = $self->{compiler};
	my $machine = $self->{machine};
	my $mpilib = $self->{mpilib};
	my $debug = $self->{debug};

	my $cmroot = $self->{configmachinesroot};
	#my @modulenodes = $cmroot->findnodes("//machine[\@MACH=\'$machine\']/module_system/modules[not(\@\*) or \@compiler=\'$compiler\' or \@mpilib=\'$mpilib\' or \@debug=\'$debug\']");
    #my @modulenodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system/modules[\@compiler=\'$compiler\' or \@mpilib=\'$mpilib\' or \@debug=\'$debug\']/module");
    #my @modulenodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system/modules[\@compiler=\'$compiler\']");
	#print Dumper \@modulenodes;
	my @modulenodes; 
	my $seqnum = 1;
	my @allmodulenodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system/modules");
	foreach my $mod(@allmodulenodes)
	{
		# If a module set does NOT have attributes, then the intention is that 
		# these modules will be loaded regardless of compiler, mpilib, or debug status.  
		if(!$mod->hasAttributes())
		{
			my @modchildren = $mod->getChildNodes();
			foreach my $child(@modchildren)
			{
				my $action = $child->getName();	
				my $actupon = $child->textContent();
				my $modhash = {action => $action, actupon => $actupon, seqnum => $seqnum};
				push(@modulenodes, $modhash);
				$seqnum += 1;
			}
		}
		else
		{
			my $attrmatch = 0;
			foreach my $qualifier (qw/ machine compiler mpilib debug /)
			{
				if($mod->hasAttribute($qualifier) && $mod->getAttribute($qualifier) eq $self->{$qualifier})
				{
				    #print "qualifier match : $qualifier\n";
					#print "mod attribute: ", $mod->getAttribute($qualifier), "\n";
					#print "self qualifier: ", $self->{$qualifier}, "\n";
					$attrmatch = 1;
				}
				elsif( $mod->hasAttribute($qualifier) && $mod->getAttribute($qualifier) ne $self->{$qualifier})
				{
				    #print "qualifier no match: $qualifier\n";
					#print "mod attribute: ", $mod->getAttribute($qualifier), "\n";
					#print "self qualifier: ", $self->{$qualifier}, "\n";
					$attrmatch = 0;
					last;
				}
			}
			if($attrmatch == 1)
			{
				my @modchildren = $mod->getChildNodes();
				foreach my $child(@modchildren)
				{
					my $action = $child->getName();
					my $actupon = $child->textContent();	
					my $modhash = { action => $action, actupon => $actupon, seqnum => $seqnum };
					push(@modulenodes, $modhash);
					$seqnum += 1;
				}
			}
		}
	}
	$self->{modulestoload} = \@modulenodes;
	#print Dumper $self;
	return @modulenodes;
}

sub writeXMLFileForCase()
{
	my $self = shift;
	#$self->moduleInit();
	my $cmroot = $self->{configmachinesroot};
	my $machine = $self->{machine};
    my @modulenodes;
    my $seqnum = 1;
    #my @allmodulenodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system");	
    my @machinenodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']");	
	
	my $machinenode = $machinenodes[0];
	
	my $casexml = XML::LibXML::Document->new("1.0.0");
	print "machinenode\n";
	print Dumper $machinenode;
	my $newmachnode = XML::LibXML::Element->new($machinenode->nodeName);
	$newmachnode->setAttribute("MACH", $machinenode->getAttribute("MACH"));
	
	my @modulesystemnodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system");
	foreach my $mnode(@modulesystemnodes)
	{
		$newmachnode->addChild($mnode);
	}

	my @envnodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']/environment_variables");
	foreach my $enode(@envnodes)
	{
		$newmachnode->addChild($enode);
	}
	my @limitnodes = $cmroot->findnodes("/config_machines/machine[\@MACH=\'$machine\']/limits");
	foreach my $lnode(@limitnodes)
	{
		$newmachnode->addChild($lnode);
	}

	my $newdom = XML::LibXML::Document->new("1.0.0");
	$newdom->setDocumentElement($newmachnode);
	my $filepath = $self->{caseroot} . "/mach_specific.xml";
	$newdom->toFile($filepath, 1);
}

sub loadModulesForCase()
{
	my $self = shift;
	my $compiler = $self->{compiler};
	my $machine = $self->{machine};
	my $mpilib = $self->{mpilib};
	my $debug = $self->{debug};

	my $cmroot = $self->{configmachinesroot};
	my @allmodulenodes = $cmroot->findnodes("/machine[\@MACH=\'$machine\']/module_system/modules");

	my @foundmodules = $self->findModules(\@allmodulenodes);
}

sub findModules()
{
	my $self = shift;
	my $modulenodes = shift;
	
	my $seqnum = 1;
	my @foundmodules;
	foreach my $mod(@$modulenodes)
	{
		# If the <modules> element doesn't have any attributes, 
		# then we want to load these modules no matter what.  
		if(!$mod->hasAttributes())
		{
			my @modchildren = $mod->getChildNodes();
			
			# for ever child node we find, 
			# action is the module action to take, actupon is the module 
			# we want to act upon, and the seqnum denotes the order in which the
			# module will be loaded. 
			foreach my $child(@modchildren)
			{
				my $action = $child->getName();
				my $actupon = $child->textContent();
				my $modhash = { action => $action, actupon => $actupon, seqnum => $seqnum} ;
				push(@foundmodules, $modhash);
				$seqnum += 1;
			} 
		}
		else
		{
			my $attrmatch = 0;
			# If the modules element has attributes, we only want to load the modules within 
			# if the attribute exists as an attribute in the modules tag, and if the attribute matches
			# the machine, compiler, etc. 
			foreach my $qualifier ( qw/ machine compiler mpilib debug/)
			{
				
				if($mod->hasAttribute($qualifier) && $mod->getAttribute($qualifier) eq $self->{$qualifier})
				{
					$attrmatch = 1;
				}
				# if the qualifier exists as an attribute but doesn't match, skip the entire block. 
				elsif( $mod->hasAttribute($qualifier) && $mod->getAttribute($qualifier) ne $self->{$qualifier})
				{
					$attrmatch = 0;
					last;
				}
			}
			if($attrmatch == 1)
			{
				my @modchildren = $mod->getChildNodes();
				foreach my $child(@modchildren)
				{
					my $action = $child->getName();
					my $actupon = $child->textContent();	
					my $modhash = { action => $action, actupon => $actupon, seqnum => $seqnum };
					push(@foundmodules, $modhash);
				}
			}
		}
	}
	return @foundmodules;
}
1;
