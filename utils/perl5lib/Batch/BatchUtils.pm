#!/usr/bin/env perl;
#==============================================================================
# File:  BatchUtils.pm
# Purpose:  Utility class for submitting jobs, and managing dependencies.  
#
#==============================================================================
use strict;
use warnings;
package Batch::BatchUtils;
use Data::Dumper;
use Cwd;
use Exporter qw(import);
use XML::LibXML;
require Batch::BatchMaker;
#push(@INC, "/gpfs/mira-fs0/projects/CESM_Atmos/usr/shollenb/devsandboxes/cesm1_3_workflow_batch/cime/utils/perl5lib/");
use lib '.';
#==============================================================================
# Base class constructor.  required args are the case name, caseroot, cime/cesmroot, 
# compiler, machine, machine root directory, the mpi library, and the scriptsroot.  
# get the paths to the config_machines and config_batch xml files, and figure out 
# the batch system type. 
#==============================================================================
sub new
{
	my ($class, %params) = @_;
	my $self = {
		case => $params{'case'}     || undef,
		caseconfig => $params{'caseconfig'} || undef,
		caseroot => $params{'caseroot'}     || undef,
		ccsmroot => $params{'ccsmroot'}     || undef,
		compiler => $params{'compiler'}  || undef,
		machine => $params{'machine'}  || undef,
		machroot => $params{'machroot'}  || undef,
		mpilib => $params{'mpilib'}  || undef,
		scriptsroot => $params{'scriptsroot'}  || undef,
	};
	bless $self, $class;

    my $perl5libdir = undef;
	my $configbatch = $self->{'machroot'} . "/config_batch.xml";
	$self->{'configbatch'} = $configbatch;
	my $configmachines = $self->{'machroot'} . "/config_machines.xml";
	$self->{'configmachines'} = $configmachines;
    my $casetoolsdir = $self->{'caseroot'} . "/Tools";
    push(@INC, $casetoolsdir);
    my $xml = XML::LibXML->new(no_blanks => 1);
    my $machineconfig = $xml->parse_file($configmachines);
    my $root = $machineconfig->getDocumentElement();
	
    my @batchtypes = $root->findnodes("/config_machines/machine[\@MACH=\'$self->{machine}\']/batch_system");
    if(! @batchtypes)
    {
        die "Could not determine batch system type for machine $self->{machine}";
    }
    $self->{'batchtype'} = $batchtypes[0]->getAttribute('name');

	$self->{dependencyqueue} = undef;
	return $self;
}
sub _check()
{
    my $self = shift;
    return 1;
}
#==============================================================================
# Get the batch system type for this machine. 
#==============================================================================
sub getBatchSystemType()
{
	my $self = shift;
	my $configmachines = $self->{'machroot'} . "/config_batch.xml";
	my $casetoolsdir = $self->{'caseroot'} . "/Tools";
	push(@INC, $casetoolsdir);
	my $xml = XML::LibXML->new(no_blanks => 1);
	my $machineconfig = $xml->parse_file($configmachines);
	my $root = $machineconfig->getDocumentElement();
	my @batchtypes = $root->findnodes("/config_machines/machine[\@MACH=\'$self->{machine}\']/batch_system");
	if(! @batchtypes)
	{
		die "Could not determine batch system type for machine $self->{machine}";
	}
	$self->{'batchtype'} = $batchtypes[0]->getAttribute('name');
}
#==============================================================================
# Get the depend_string so jobs can be submitted with dependencies.
#==============================================================================
sub getDependString()
{
	my $self = shift;
	my $jobid = shift;
	my $xml = XML::LibXML->new(no_blanks => 1);
	my $batchconfig = $xml->parse_file($self->{'configbatch'});
    my $root = $batchconfig->getDocumentElement();
	my @dependargs = $root->findnodes("/config_batch/batch_system[\@type=\'$self->{'batchtype'}\']/depend_string");
	if(! @dependargs)
	{
		die "could not find depend string for this batch system type";
	}
	my $deparg = $dependargs[0]->textContent();
    $deparg =~ s/jobid/$jobid/g;
	return $deparg;
}

#==============================================================================
# Get the job id from the output of the job submission. 
#==============================================================================
sub getJobID()
{
	my $self = shift;
	my $jobstring = shift;
	chomp $jobstring;
	my $xml = XML::LibXML->new(no_blanks => 1);
    my $batchconfig = $xml->parse_file($self->{'configbatch'});
    my $root = $batchconfig->getDocumentElement();
	my @jobidpatterns = $root->findnodes("/config_batch/batch_system[\@type=\'$self->{'batchtype'}\']/jobid_pattern");
	if(!@jobidpatterns)
	{
		die "could not find job id pattern for batch system type $self->{'batchtype'}";
	}
	
	my $jobid = undef;
	my $jobidpat = $jobidpatterns[0]->textContent();
	my $pattern  = qr($jobidpat);
	if($jobstring =~ /$pattern/ )
	{
		$jobid = $1;
	}
	return $jobid;
}


#==============================================================================
# Submit jobs in a dependency chain.
#==============================================================================
sub submitJobs()
{
	my $self = shift;
	my $depjobid = undef;
	
	my %depqueue = %{$self->{dependencyqueue}};
	my $lastjobseqnum = (sort {$b <=> $a } keys %depqueue)[0];
	foreach my $jobnum(sort keys %depqueue)
	{
		foreach my $jobname(@{$depqueue{$jobnum}})
		{
			my $islastjob = 0;
			$islastjob = 1 if ($jobnum == $lastjobseqnum);
			$depjobid = $self->submitSingleJob($jobname, $depjobid, $islastjob);
		}
	}
}

#==============================================================================
# Submit a job, given the script name. If a dependent job id is given, 
# use it as the dependency. If the islastjob flag is passed in, add the 
# islastjob=TRUE flag to the environment of the job submission.. 
#==============================================================================
sub submitSingleJob()
{
	my $self = shift;
	my $scriptname = shift;
	my $dependentJobId = shift;
	my $islastjob = shift;
	my %config = %{$self->{'caseconfig'}};
	my $dependarg = '';
	my $submitargs = '';
    $submitargs = $self->getSubmitArguments($scriptname, $dependentJobId);
	if(! defined $submitargs && length($submitargs) <= 0)
	{
	    $submitargs = '' ;
	}
	#if(defined $dependentJobId)
	#{
	#	#$dependarg = $self->getDependString($dependentJobId);	
    #    $submitargs = $self->getSubmitArguments($dependentJobId);
	#	#print "depend arg: $dependarg\n";
    #    print "submit args: $submitargs\n";
	#}
	if($islastjob)
	{
		$ENV{'islastjob'} = 'TRUE';
	}
	print "Submitting CESM job script $scriptname\n";
	my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname";
    
	print "submitting command $runcmd\n";
	open (my $RUN, "-|", $runcmd) or die "cannot run the command \"$runcmd\"";

	my $output = <$RUN>;
	chomp $output;	
	print "submit output: |$output|\n";
	
    #print "last dependent job number: $lastjobseqnum\n";
    #print "last dependent job: ", @{$depqueue{$lastjobseqnum}} ,  "\n";
	
	my $jobid = $self->getJobID($output);
	return $jobid;
}


#==============================================================================
# Check to see if the next set of jobs needs to be submitted.  
#==============================================================================
sub doResubmit()
{
	my $self = shift;
    my $islastjob = shift;
    my $resubmit = shift;
	
    # If the islastjob flag is true, and resubmit is > 0,  do the dependency
    # check and job resubmission again 
    if($islastjob eq 'TRUE' && $resubmit > 0)
    {
	    my %config = %{$self->{caseconfig}};
	    $self->dependencyCheck();
	    $self->submitJobs();		
	    my $newresubmit = $config{'RESUBMIT'} - 1;
	    my $owd = getcwd;
	    chdir $config{'CASEROOT'};
	    `./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit`;
	    if($?)
	    {
	    	print "could not execute ./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit\n";
	    }
	    chdir $owd;
    }
    else    
    {
        return;
    }
	
}

#==============================================================================
# If we need to resubmit a CESM job + post-run jobs, this subroutine will check the 
# env*.xml variables to see which jobs need to be resubmitted.  
#==============================================================================
sub dependencyCheck()
{
	my $self = shift;
	my %config = %{$self->{'caseconfig'}};
	# we always want to run the CESM test or run again..
	if(-e "$config{'CASE'}.test")
	{
		my $jobname = "$config{'CASE'}.test";
		$self->addDependentJob($jobname);
        return;
	}
	else
	{
		my $jobname = "$config{'CASE'}.run";
		$self->addDependentJob($jobname);
	}
	
	# do we add the short-term archiver to the dependency queue? 
	if($config{'DOUT_S'} eq 'TRUE')
	{
		my $jobname = "$config{'CASE'}.st_archive";
		$self->addDependentJob($jobname);
	}
	
}

#==============================================================================
# Adds a job script name or array of job script names to the dependent jobs queue. 
#==============================================================================
sub addDependentJob()
{
	my $self = shift;
	# Either a string with the job name or an array of job names
	my $jobref = shift;
	my $jobcounter = 0;
	
	# set up the dependency hash if not done.  
	if(! defined $self->{dependencyqueue})
	{
		$self->{dependencyqueue} = {};
	}
	# get the dependency hash. 
	my %dependencyqueue = %{$self->{dependencyqueue}};
	
	# Increment the job counter for each job set in the dependency queue. 
	foreach my $jobnum(keys %{$self->{dependencyqueue}})
	{
		$jobcounter += 1;
	}

	# If the jobref is a regular string scalar, make an array out of it 
	# and add it to the dependency queue. 
	if(! ref($jobref) )
	{
		my @jobarray = ( $jobref );
		$dependencyqueue{$jobcounter} = \@jobarray;
	}
	# if we have an array, add the array to the dependency queue directly.  
	elsif(ref($jobref) eq 'ARRAY')
	{
		$dependencyqueue{$jobcounter} = $jobref;
	}
	$self->{dependencyqueue} = \%dependencyqueue;
}

#==============================================================================
#==============================================================================
sub getSubmitArguments()
{
    my $self = shift;
    # We need the script name and the dependent job id.
    my $scriptname = shift;
    my $dependentjobid = shift;

    my $batchmaker = Batch::BatchFactory::getBatchMaker( caseroot => $self->{caseroot}, case => $self->{case},
                                                  mpilib => $self->{mpilib}, scriptsroot => $self->{scriptsroot},
                                                  machroot => $self->{machroot}, machine => $self->{machine},
                                                  compiler => $self->{compiler} );


    # Find the submit arguments for this particular batch system.
    my $xml = XML::LibXML->new(no_blanks => 1);
    my $batchconfig = $xml->parse_file($self->{'configbatch'});
    my $root = $batchconfig->getDocumentElement();

    my @dependargs = $root->findnodes("/config_batch/batch_system[\@type=\'$self->{'batchtype'}\']/submit_args/arg");

    my $submitargs = '';

    if(@dependargs)
	{
		foreach my $dependarg(@dependargs)
    	{

    	    my $argFlag = $dependarg->getAttribute('flag');
    	    my $argName = $dependarg->getAttribute('name');
    	    if(defined $argName && length($argName) > 0)
    	    {
    	        # Get the actual data field from the BatchMaker class.
    	        my $field = $batchmaker->getField($argName);
    	        if(! defined $field)
    	        {
    	            die "$argName not defined! Aborting...";
    	        }
    	        else
    	        {
    	            $submitargs .= " $argFlag $field";
    	        }
    	    }
    	    # If the argName isn't defined,
    	    elsif(defined $argFlag && ! defined $argName)
    	    {
    	        $submitargs .= " $argFlag";
    	    }
    	}
	}

    if(defined $dependentjobid)
	{
		my $dependArg = $self->getDependString($dependentjobid);
		$submitargs .= " $dependArg ";
	}

	return $submitargs;

}


#==============================================================================
# Factory package for getting the correct BatchUtils class.  Either machine-specific or 
# batch system specific classes can be defined in this package, and that specific
# class will be returned if it exists, or the base class will be returned if 
# nothing machine-specific or batch system specific exists.  
#==============================================================================
package Batch::BatchUtilsFactory;
use Exporter qw(import);
use XML::LibXML;
use Data::Dumper;
sub getBatchUtils
{
    my (%params) = @_;
    
    # We need a machine to be defined
    my $machine = $params{'machine'};
	if(!defined $machine)
	{
		die "BatchUtilsFactory: machine must be defined!";
	}
    # Find the batch system type based on the machine. 
	my $batchtype = getBatchSystemType($params{'machine'}, $params{'machroot'}, $params{'caseroot'});

    # Make a new base class 
    my $batchutils = Batch::BatchUtils->new(%params);
    
    # Get the machine-specific and batch-specific class names
    my $machclassname = "Batch::BatchUtils_" . $machine; 
    my $batchclassname = "Batch::BatchUtils_" . $batchtype;

    # Try to re-bless the base class into the machine-specific class. 
    # then test to see if we can actually call methods. If the 
    # eval fails, then we know the machine-specific class doesn't exist. 
    my $rv = eval 
    { 
        #$machclassname->new(%params);
        #require $machclassname;
        bless $batchutils, $machclassname;
        $batchutils->_test();
        1;
    };
    # If we don't get an error, return the machine-specific 
    # BatchUtils object. 
    if(! $@)
    {
        return $batchutils;
    }
    else
	{
		bless $batchutils, "Batch::BatchUtils";
	}

    # Now try to create the batch-system specific class. 
    $rv = eval 
    { 
        #$batchclassname->new(%params);
        #require $batchclassname;
        bless $batchutils, $batchclassname; 
        $batchutils->_test();
        #1;
    };

    # Return the batch-system specific class if we don't get an error
    if(! $@)
    {
        return $batchutils;
    }
    else    
    {
		#print "batch-specific error found, $@";
		bless $batchutils, "Batch::BatchUtils";
        return $batchutils;
    }

}

#==============================================================================
# BatchUtilsFactory getBatchSystemType method. 
#==============================================================================
sub getBatchSystemType()
{
    my $machine = shift;
    my $machroot = shift;
    my $caseroot = shift;
    my $configbatch = "$machroot/config_batch.xml";
    my $configmachines = "$machroot/config_machines.xml";
    my $casetoolsdir = "$caseroot/Tools";
    push(@INC, $casetoolsdir);
    my $xml = XML::LibXML->new(no_blanks => 1);
    my $batchconfig = $xml->parse_file($configbatch);
    my $machineconfig = $xml->parse_file($configmachines);
    my $root = $machineconfig->getDocumentElement();
    my @batchsystems = $root->findnodes("/config_machines/machine[\@MACH=\'$machine\']/batch_system");
    if(! @batchsystems)
    {
        die "Could not determine batch system type for machine $machine";
    }
    #my $batchtype = $batchsystems[0]->textContent();
    my $batchtype = $batchsystems[0]->getAttribute('name');
    #print "Found batch type $batchtype\n";
    return $batchtype;
}


#==============================================================================
# Mira/ALCF specific BatchUtils class, since the workflow for ALCF has to be 
# completely different.  
#==============================================================================
package Batch::BatchUtils_mira;
use base qw( Batch::BatchUtils );
use Data::Dumper;
use Cwd;

#==============================================================================
# Overridden submitJobs() method for Mira.  
#==============================================================================
sub submitJobs()
{
    my $self = shift;
    my $depjobid = shift;
    my %depqueue = %{$self->{dependencyqueue}};

    # Get the first job sequence number. 
    my $firstjobseqnum = (sort {$a <=> $b } keys %depqueue)[0];
    # we get the first job array reference. 
    my $firstjobarray = $depqueue{$firstjobseqnum};
    # Get the first job name. 
    my $firstjobname = $$firstjobarray[0];

    # submit the CESM run, and nothing else. 
    $depjobid = $self->submitSingleJob($firstjobname, $depjobid, 0);
}

sub submitSingleJob()
{
    my $self = shift;
    my $scriptname = shift;
    my $dependentJobId = shift;
    my $islastjob = shift;
    my $workflowhostfile = "./workflowhostfile";
    if(! -e $workflowhostfile)
    {
        open (my $W, ">", $workflowhostfile) or die "could not open workflow host file, $!";
        if(defined $ENV{'HOST'})
        {
            print $W $ENV{'HOST'} . "\n";
        }
        elsif(defined $ENV{'HOSTNAME'})
        {
            print $W $ENV{'HOSTNAME'} . "\n";
        }
        close $W;
    }
    $self->SUPER::submitSingleJob($scriptname, $dependentJobId, $islastjob);

}
#==============================================================================
# Mira-specific doResubmit call.  If this is called from the cesm run, then we 
# have to ssh over to tukey and run the short-term archiver. 
#==============================================================================
sub doResubmit()
{
    my $self = shift;
    my $islastjob = shift;
    my $resubmit = shift;
    my $scriptname = shift;

    print "in doResubmit: scriptname: $scriptname\n";
    my %config = %{$self->{'caseconfig'}};
    if($scriptname =~ /run/)
    {
        my $starchivescript = $scriptname;
        $starchivescript =~ s/run/st_archive/g;
        
        my $submitargs = $self->getSubmitArguments($starchivescript);
        if($config{'RESUBMIT'} > 0 && $config{'CONTINUE_RUN'} eq 'TRUE')
        {
            $submitargs .= " --env islastjob=TRUE ";
        }
            
        my $submitstuff = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $starchivescript";
        #my $cmd = "ssh tukeylogin1 qsub  -A CESM_Atmos -t 60 -n 1 -q default --mode script ./$starchivescript";
        my $runcmd = "ssh tukeylogin1 $submitstuff";
        print "runcmd is: $runcmd\n";
        qx($runcmd) or die "could not exec cmd $runcmd, $!";
    }

    if($scriptname =~ /archive/ && $islastjob eq 'TRUE' && $resubmit > 0)
    {
        my $newresubmit = $config{'RESUBMIT'} - 1;
        my $owd = getcwd;
        chdir $config{'CASEROOT'};
        `./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit`;
        chdir $owd;
        
        my $runscript = $scriptname;
        $runscript =~ s/st_archive/run/g;
        
        my $submitargs = $self->getSubmitArguments($runscript);
        my $submitstuff = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $runscript";
        open (my $W, "<", "./workflowhostfile" ) or die "could not open workflow host file, $!";
        my $text = <$W>;
        close $W;
        my $runhost;
        if($text =~ /mira/)
        {
            $runhost = "miralac1";
        }
        elsif($text =~ /cetus/)
        {
            $runhost = "cetuslac1";
        }
        my $runcmd = "ssh tukeylogin1 ssh $runhost $submitstuff ";
        qx($runcmd) or die "could not exec cmd $runcmd, $!";
        
    }
    
}

# If we need arguments when we submit a job, this is where we will add them.
# This is really specific to Mira/ALFC, no center requires qsub submit arguments. 
sub getSubmitArguments()
{
    my $self = shift;
    
    # We need the script name and the dependent job id. 
    my $scriptname = shift;
    my $dependentjobid = shift;

    my $batchmaker = Batch::BatchFactory::getBatchMaker( caseroot => $self->{caseroot}, case => $self->{case},
                                                  mpilib => $self->{mpilib}, scriptsroot => $self->{scriptsroot},
                                                  machroot => $self->{machroot}, machine => $self->{machine},
                                                  compiler => $self->{compiler} );
    
    # Set the node count to 1 if this is the short-term archive script. 
    if(defined $scriptname && $scriptname =~ /archive/)
    {
        $batchmaker->overrideNodeCount(1);
    }

    # Find the submit arguments for this particular batch system. 
    my $xml = XML::LibXML->new(no_blanks => 1);
    my $batchconfig = $xml->parse_file($self->{'configbatch'});
    my $root = $batchconfig->getDocumentElement();

    my @dependargs = $root->findnodes("/config_batch/batch_system[\@type=\'$self->{'batchtype'}\']/submit_args/arg");

    my $submitargs = '';

    # get the flag, and the name of the flag.  
    foreach my $dependarg(@dependargs)
    {

        my $argFlag = $dependarg->getAttribute('flag');
        my $argName = $dependarg->getAttribute('name');
        #print "flag: $argFlag\n";
        if(defined $argName && length($argName) > 0)
        {
            # Get the actual data field from the BatchMaker class. 
            my $field = $batchmaker->getField($argName);
            if(! defined $field)
            {
                die "$argName not defined! Aborting...";
            }
            else
            {
                $submitargs .= " $argFlag $field";
            }
        }
        # If the argName isn't defined, 
        elsif(defined $argFlag && ! defined $argName)
        {
            $submitargs .= " $argFlag";
        }
    }

    if(defined $dependentjobid)
    {
        my $dependArg = $self->getDependString($dependentjobid);
        $submitargs .= " $dependArg";
    }
    # Need to add the --cwd argument to the submit args if this is the st_archive script
    if(defined $scriptname && $scriptname =~ /archive/)
    {
        $submitargs .= " --cwd $self->{'caseroot'} ";
    }
        
    #print "submit args are now: $submitargs\n";
    return $submitargs;
}

sub _test()
{
    my $self = shift;
}

1;
