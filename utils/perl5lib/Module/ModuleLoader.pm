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

	#ndk: to help with debugging, I was creating a timestamped file
        #     everytime this function was called.  Set $debugML=1 to turn on.
        my $debugML=0;
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
        my ($ntimestamp, $dfile);
	if ($debugML) {
           ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
           $ntimestamp = sprintf ( "%04d.%02d.%02d.%02d.%02d.%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
           $dfile="debugML_".$ntimestamp.".txt";
           open(DEBUGOUT, ">$dfile");
	}

	# Disable check - check will cause problem if scheduler is changing env - better reload 
	#if(defined $ENV{CIME_MODULES_LOADED}) {
	#    if($debugML) print DEBUGOUT " ENV CIME_MODULES_LOADED is true, returning\n";
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
        if ($debugML) {
           print DEBUGOUT "time=$ntimestamp running command $cmd\n";
        }

        my @output;
	eval { @output = qx($cmd);};
	chomp @output;

        if ($debugML) {
           print DEBUGOUT "after running command output follows\n";
           print DEBUGOUT "number of lines = $#output\n";
           print DEBUGOUT '@output=qx($cmd) and output='; 
           print DEBUGOUT "\n\n";
           for my $i(0 .. $#output) {
              my $oo=$output[$i];
              print DEBUGOUT "$oo\n";
           }
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
                    if ($debugML) {print DEBUGOUT "line: $line\n";}
	            my ($key, $value) = split('=', $line, 2);
	            $newenv{$key} = $value;
                }
	}
	my %newbuildenv;

	# ndk: add a loop to see if our above module adjusting _removed_ 
        # any environment variables which may need to be accounted for
	foreach my $k(keys %oldenv) {
           if($k =~ /BASH_FUNC/i) {
              # leave this particular one alone
           } elsif (defined $oldenv{$k} && !defined $newenv{$k}) {
	      # if key in old but NOT in new, consider removing
              if ($debugML) {print DEBUGOUT "del $k=$oldenv{$k}\n";}
	      delete $ENV{$k};
           }
        }
        
        foreach my $k(keys %newenv)
	{
           if ($debugML) {print DEBUGOUT "k=$k ";}
           if(! defined $oldenv{$k})   # if this key is _not_ in the old set, add it as new
           {
              if ($debugML) {print DEBUGOUT "new $k=$newenv{$k}\n";}
              $newbuildenv{$k} = $newenv{$k};
              $ENV{$k} = $newenv{$k};
           }
           elsif(defined $oldenv{$k} && $newenv{$k} ne $oldenv{$k}) 
           {
              # if this var is in the old, and is different than it was in old
              if ($debugML) {print DEBUGOUT "mod $k=$newenv{$k}\n";}
              $newbuildenv{$k} = $newenv{$k};
              $ENV{$k} = $newenv{$k};
           }
           else {
              if ($debugML) {print DEBUGOUT "\n";}
           }
	}

        if ($debugML) {
           foreach my $k(keys %newbuildenv) {
              print DEBUGOUT "nbe $k=$newbuildenv{$k}\n";
           }
           foreach my $k(keys %ENV) {
              print DEBUGOUT "ENV $k=$ENV{$k}\n";
           }
           close (DEBUGOUT);
        }

	$ENV{CIME_MODULES_LOADED} = 1;
	return %newbuildenv;

}
1;
