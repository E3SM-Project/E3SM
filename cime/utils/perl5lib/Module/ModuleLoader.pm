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

        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
        my $ntimestamp = sprintf ( "%04d.%02d.%02d.%02d.%02d.%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
        #print "ntimestamp=$ntimestamp\n";
        my $dfile="debugML_".$ntimestamp.".txt";
        open(DEBUGOUT, ">$dfile");

	# Disable check - check will cause problem if scheduler is changing env - better reload 
	#ndk if(defined $ENV{CIME_MODULES_LOADED}) {return $self};
	#if(defined $ENV{CIME_MODULES_LOADED}) {
	#    print "ndk neon env CIME_MODULES_LOADED is true, returning\n";
	#    print DEBUGOUT "ndk neon env CIME_MODULES_LOADED is true, returning\n";
	#    return $self
	#};

	my %oldenv = %ENV;
	my $envfile = $self->{'caseroot'} . "/env_mach_specific";
	my $cshenv = "env ";
    	$cshenv .= "COMPILER=". $self->{'compiler'} . " " ;
    	$cshenv .= "MPILIB=". $self->{'mpilib'} . " " ;
    	$cshenv .= "DEBUG=". $self->{'debug'} . " " ;
    	$cshenv .= "CASEROOT=". $self->{'caseroot'} . " " ;
    	$cshenv .= "PERL=TRUE ";
	my $cmd = $cshenv . " " . $envfile ;
	#print "running command $cmd\n";
	print "ndk neon running command $cmd\n";
        print DEBUGOUT "ndk neon time=$ntimestamp running command $cmd\n";
	my @output;
	eval { @output = qx($cmd);};  #ndk looking at this while debugging problems with modules on NERSC systems with bash
	#eval { @output = qx(env $envfile);};
        #eval { @output = qx(env COMPILER=intel MPILIB=mpt DEBUG=FALSE CASEROOT=/scratch2/scratchdirs/ndk/acme_scratch/SMS.ne30_oEC.A_WCYCL2000.edison_intel.b07-feb25 PERL=TRUE $envfile);};
        #eval { @output = qx(env COMPILER=intel MPILIB=mpt DEBUG=FALSE PERL=TRUE $envfile);};
        #eval { @output = qx(env COMPILER=intel MPILIB=mpt DEBUG=FALSE $envfile);}; #works
        #eval { @output = qx(env COMPILER=intel MPILIB=mpt DEBUG=FALSE CASEROOT=/scratch2/scratchdirs/ndk/acme_scratch/SMS.ne30_oEC.A_WCYCL2000.edison_intel.b07-feb25 $envfile);}; #works
	#eval { $out = `$cmd`;};
	#$out = `$cmd`;
	chomp @output;
        print DEBUGOUT "after running command output follows\n";
        print DEBUGOUT "number of lines = $#output\n";
        print DEBUGOUT '@output=qx($cmd) and output='; 
        print DEBUGOUT "\n\n";
        for my $i(0 .. $#output) {
           my $oo=$output[$i];
           print DEBUGOUT "$oo\n";
        }
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
                    print DEBUGOUT "line: $line\n";
	            my ($key, $value) = split('=', $line, 2);
	            $newenv{$key} = $value;
                }
	}
	my %newbuildenv;

        # only thing i'm suggesting we add in this case -- the other edits are debug
	if (defined $oldenv{'_LMFILES_'}  && defined $newenv{'_LMFILES_000'} ) { # if there is a LF in old and a LF000 in new
            # _LMFILES_ must have reached max size and was split into _LMFILES_000 and _LMFILES_001, ...
	    delete $ENV{'_LMFILES_'}
	}

	foreach my $k(keys %newenv)
	{
	    print DEBUGOUT "k=$k ";
		if(! defined $oldenv{$k})   # if this key is _not_ in the old set, add it as new
		{
		    #print "new env var: $newenv{$k}\n";
                   print DEBUGOUT "new $k=$newenv{$k}\n";
			$newbuildenv{$k} = $newenv{$k};
			$ENV{$k} = $newenv{$k};
		}
		elsif(defined $oldenv{$k} && $newenv{$k} ne $oldenv{$k}) # if this var is in the old, and is different than it was in old
		{
		    #print "modified env var: $newenv{$k}\n";
                    print DEBUGOUT "mod $k=$newenv{$k}\n";
			$newbuildenv{$k} = $newenv{$k};
			$ENV{$k} = $newenv{$k};
		}
	    else{
		print DEBUGOUT "\n";
	    }
	}

	foreach my $k(keys %newbuildenv) {
	    print DEBUGOUT "nbe $k=$newbuildenv{$k}\n";
	}
	foreach my $k(keys %ENV) {
	    print DEBUGOUT "ENV $k=$ENV{$k}\n";
	}
        close (DEBUGOUT);

	$ENV{CIME_MODULES_LOADED} = 1;
	return %newbuildenv;

}
1;
