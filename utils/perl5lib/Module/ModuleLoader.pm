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
	
	my @initnodes = $xml->findnodes("//machine[\@MACH=\'$machine\']/module_system/init_path[\@lang=\'perl\']");
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

	#print "new mach node: ", $newmachnode->toString();
	my $newdom = XML::LibXML::Document->new("1.0");
	$newdom->setDocumentElement($newmachnode);
	my $filepath = $self->{caseroot} . "/mach_specific.xml";
	$newdom->toFile($filepath, 1) || die "could not write file: ", $self->{caseroot} . ", $?\n";
}

sub findModulesForCase()
{
	my $self = shift;
	my $compiler = $self->{compiler};
	my $machine = $self->{machine};
	my $mpilib = $self->{mpilib};
	my $debug = $self->{debug};

	#my $cmroot = $self->{configmachinesroot};
	my $machspecificfile = $self->{'caseroot'} . "/mach_specific.xml";
	if( ! -e $machspecificfile)
	{
		die "$machspecificfile was not found!\n";
	}
    my $parser = XML::LibXML->new(no_blanks => 1);
    #my $xml = $parser->parse_file($self->{'configmachinesfile'});
	my $casemoduleparser = $parser->parse_file($machspecificfile);
	my @allmodulenodes = $casemoduleparser->findnodes("/machine[\@MACH=\'$machine\']/module_system/modules");

	my @foundmodules = $self->findModules(\@allmodulenodes);
	$self->{modulestoload} = \@foundmodules;
	#print Dumper $self;
	return @foundmodules; 
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
				    $seqnum += 1;
				}
			}
		}
	}
	return @foundmodules;
}

sub loadModules()
{
	my $self = shift;
	
	if(! defined $self->{modulestoload})
	{
		die "no modules to load, wtf??";
	}
	
	my $modulestoload = $self->{modulestoload};
	
	foreach my $mod(@$modulestoload)
	{
		print "mod seqnum: $mod->{seqnum}\n";
		print "mod action: $mod->{action}\n";
		print "mod actupon: $mod->{actupon}\n";
		my $cmd = $self->{cmdpath} . " $mod->{action}  $mod->{actupon}";
		print "running cmd: $cmd\n";
		eval qx($cmd);
		if($?)
		{
			warn "module cmd $cmd died with $? $!\n";
		}
	}
	#map { print "key: $_, value: $ENV{$_}\n" } sort keys %ENV;
	$self->writeCshModuleFile();
	$self->writeBashModuleFile();
}

sub getEnvironmentVars()
{
	my $self = shift;
    #my @envnodes = $xml->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system/environment_variables/env");
	#
	#foreach my $enode(@envnodes)
	#{
	#	my $type = $enode->getAttribute('name');
	#	my $value = $enode->getValue();
	#}
}

sub getLimits()
{
	my $self = shift;
}

sub writeCshModuleFile()
{
	my $self = shift;
	my $machine = $self->{machine};
	
	#my $xml = $self->{configmachinesroot};
    my $parser = XML::LibXML->new(no_blanks => 1);
    my $xml = $parser->parse_file($self->{'configmachinesfile'});
	
       #my @initnodes = $xml->findnodes("//machine[\@MACH=\'$machine\']/module_system/init_path[\@lang=\'perl\']");	
	my @cshinitnodes = $xml->findnodes("//machine[\@MACH=\'$machine\']/module_system/init_path[\@lang=\'csh\']");
	print Dumper $xml;
	print Dumper \@cshinitnodes;
	
	die "no csh init path defined for this machine!" if !@cshinitnodes;
	foreach my $node(@cshinitnodes)
	{
		$self->{cshinitpath} = $node->textContent();
	}
	
	

my $csh =<<"START";
#!/usr/bin/env csh -f 
#===============================================================================
# Automatically generated module settings for $self->{machine}
#===============================================================================

source $self->{cshinitpath}
START
	
	my @allmodules = $xml->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system/modules");
	
	foreach my $mod(@allmodules)
	{
		if(!$mod->hasAttributes())
		{
			my @modchildren = $mod->getChildNodes();
			foreach my $child(@modchildren)
			{
				my $action = $child->getName();
				my $actupon = $child->textContent();
				$csh .= "module $action $actupon\n";
			}
		}
		else
		{
			my @attrs = $mod->attributes;
	
			$csh .= "if ( ";
			while(@attrs)
			{
				my $attr = shift @attrs;
				my $name = uc($attr->getName());
				my $value = $attr->getValue();
				$csh .= "\$$name == \"$value\"";
				$csh .= " && " if(@attrs);
			}
			$csh .= " ) then\n";

			my @modchildren = $mod->getChildNodes();
			foreach my $child(@modchildren)
			{
				
				my $action = $child->getName();
				my $actupon = $child->textContent();
				$csh .= "\tmodule $action $actupon\n";
			}
			$csh .= "endif\n";
		}
	}
	
	my @envnodes = $xml->findnodes("/config_machines/machine[\@MACH=\'$machine\']/environment_variables/env");
    foreach my $enode(@envnodes)
    {
        my $name = $enode->getAttribute('name');
        my $value = $enode->textContent();
		$csh .= "setenv $name $value\n";
    }

	my @limitnodes = $xml->findnodes("/config_machines/machine[\@MACH=\'$machine\']/limits/limit");
	foreach my $lnode(@limitnodes)
	{
		my $name = $lnode->getAttribute('name');
		my $value = $lnode->textContent();
		$csh .= "limit $name $value\n";
	}
	
	open my $CSHFILE, ">", "$self->{caseroot}/.env_mach_specific.csh" || die " coult not open test.csh, $!";
	print $CSHFILE $csh;
	close $CSHFILE;
}

sub writeBashModuleFile()
{
    my $self = shift;
	my %cshtobash = ( "cputime" => "-t", 
                      "filesize" => "-f",
                      "datasize" => "-d", 
                      "stacksize" => "-s", 
                      "coredumpsize" => "-c", 
                      "memoryuse" => "-m", 
                      "vmemoryuse" => "-v",
                      "descriptors" => "-n", 
                      "memorylocked" => "-l",
                      "maxproc" => "-u" );
    my $machine = $self->{machine};

    #my $xml = $self->{configmachinesroot};
    my $parser = XML::LibXML->new(no_blanks => 1);
    my $xml = $parser->parse_file($self->{'configmachinesfile'});

    my @bashinitnodes = $xml->findnodes("//machine[\@MACH=\'$machine\']/module_system/init_path[\@lang=\'bash\']");

    die "no bash init path defined for this machine!" if !@bashinitnodes;
    foreach my $node(@bashinitnodes)
    {
        $self->{bashinitpath} = $node->textContent();
    }

    my $bash =<<"START";
#!/usr/bin/env bash -f 
#===============================================================================
# Automatically generated module settings for $self->{machine}
#===============================================================================

.  $self->{bashinitpath}
START

	my @allmodules = $xml->findnodes("/config_machines/machine[\@MACH=\'$machine\']/module_system/modules");
	
	foreach my $mod(@allmodules)
	{
		if(! $mod->hasAttributes())
		{
			my @modchildren = $mod->getChildNodes();
			foreach my $child(@modchildren)
			{
				my $action = $child->getName();
				my $actupon = $child->textContent();
				$bash .= "module $action $actupon \n";
			}
		}
		else
		{
			my @attrs = $mod->attributes;
			
			$bash .= "if [ ";
			while(@attrs)
			{
				my $attr = shift @attrs;
				my $name = uc($attr->getName());
				my $value = $attr->getValue();
				$bash .= "\"\$$name\" = \"$value\"";
				$bash .= "] && [ " if (@attrs);
			}
			$bash .= " ]\n";
			$bash .= "then\n";
	
			my @modchildren = $mod->getChildNodes();
			foreach my $child(@modchildren)
			{
                my $action = $child->getName();
				my $actupon = $child->textContent();
				$bash .= "\tmodule $action $actupon\n";
			}
			$bash .= "fi\n";
		}
	}


    my @envnodes = $xml->findnodes("/config_machines/machine[\@MACH=\'$machine\']/environment_variables/env");
    foreach my $enode(@envnodes)
    {
        my $name = $enode->getAttribute('name');
        my $value = $enode->textContent();
        $bash .= "export $name=$value\n";
    }

    my @limitnodes = $xml->findnodes("/config_machines/machine[\@MACH=\'$machine\']/limits/limit");
    foreach my $lnode(@limitnodes)
    {
        my $name = $lnode->getAttribute('name');
        my $value = $lnode->textContent();
		print Dumper \%cshtobash;
		my $bashname = $cshtobash{$name};
        $bash .= "ulimit $bashname $value\n";
    }

	open my $BASHFILE, ">", "$self->{caseroot}/.env_mach_specific.bash" || die "could not open .env_mach_specific.bash, $!";
	print $BASHFILE $bash;
	close $BASHFILE;
}
1;
