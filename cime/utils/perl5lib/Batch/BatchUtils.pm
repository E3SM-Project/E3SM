#!/usr/bin/env perl;
#==============================================================================
# File:  BatchUtils.pm
# Purpose:  Utility class for submitting jobs, and managing dependencies.  
# PLEASE NOTE: for each subclass you define, you MUST add a _test() method 
# to the subclass for the factory to work 
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
    $self->{'batchtype'} = $batchtypes[0]->getAttribute('type');

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
	my @machjobidpatterns = $root->findnodes("/config_batch/batch_system[\@MACH=\'$self->{machine}\']/jobid_pattern");
	my $jobidpat;
	if(@machjobidpatterns)
	{
		$jobidpat = $machjobidpatterns[0]->textContent();
	}
	else
	{
		my @basejobidpatterns = $root->findnodes("/config_batch/batch_system[\@type=\'$self->{'batchtype'}\']/jobid_pattern");
		if(!@basejobidpatterns)
		{
			die "could not find job id pattern for batch system type $self->{'batchtype'}";
		}
		else
		{
			$jobidpat = $basejobidpatterns[0]->textContent();
		}
	}
	
	my $jobid = undef;
	my $pattern  = qr($jobidpat);
	if($jobstring =~ /$pattern/ )
	{
		$jobid = $1;
	}
	else
	{
		die " could not ascertain dependent job id... aborting";
	}
	return $jobid;
}


#==============================================================================
# Submit jobs in a dependency chain.
#==============================================================================
sub submitJobs()
{
	my $self = shift;
	my $sta_ok = shift;
	my $depjobid = shift;
	
	my %depqueue = %{$self->{dependencyqueue}};
	my $lastjobseqnum = (sort {$b <=> $a } keys %depqueue)[0];
	foreach my $jobnum(sort keys %depqueue)
	{
		foreach my $jobname(@{$depqueue{$jobnum}})
		{
			my $islastjob = 0;
			$islastjob = 1 if ($jobnum == $lastjobseqnum);
            $depjobid = $self->submitSingleJob($jobname, $depjobid, $islastjob, $sta_ok);
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
	my $sta_ok = shift;
	my %config = %{$self->{'caseconfig'}};
	my $dependarg = '';
	my $submitargs = '';

    $submitargs = $self->getSubmitArguments($scriptname, $dependentJobId);

	if(! defined $submitargs && length($submitargs) <= 0)
	{
	    $submitargs = '' ;
	}

	if($islastjob)
	{
		$ENV{'islastjob'} = 'TRUE';
	}
	else
	{
		$ENV{'islastjob'} = 'FALSE';
	}
	#my $sta_argument = '';
	if(defined $sta_ok)
	{
		#$sta_argument = " -F \"--sta_ok\"";
		$ENV{'sta_ok'} = 'TRUE';
	}
	else
	{
		#$ENV{'sta_ok'} = 'FALSE';
		delete $ENV{'sta_ok'};
	}
	print "Submitting CESM job script: $scriptname\n";
	#my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname $sta_argument";
	my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname ";
	print ": $runcmd\n";    
	my $output;

	eval {
        open (my $RUN, "-|", $runcmd) // die "job submission failed, $!";

		foreach (<$RUN>) {
			chomp;
			print "$_\n";
			$output .= $_;
		}

		close $RUN or die "job submission failed: |$?|, |$!|"
	};
	my $exitstatus = ($?>>8);
	if($exitstatus != 0)
	{
		print "Job submission failed\n";
		exit(1);
	}
		
	chomp $output;	
	
	my $jobid = $self->getJobID($output);
	print "Job ID: $jobid\n";
	return $jobid;
}


#==============================================================================
# Base class doResubmit
# Check to see if the next set of jobs needs to be submitted.  
#==============================================================================
sub doResubmit()
{
	my $self = shift;
    my $islastjob = shift;
    my $resubmit = shift;
	my $sta_ok = shift;
	
    # If the islastjob flag is true, and resubmit is > 0,  do the dependency
    # check and job resubmission again 
    if($islastjob eq 'TRUE' && $resubmit > 0 && defined $sta_ok)
    {
	    my %config = %{$self->{caseconfig}};
	    $self->dependencyCheck("sta_ok");
	    $self->submitJobs("sta_ok");		
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
# For now, we are only handling cesm runs and the short-term archiver. 
#==============================================================================
sub dependencyCheck()
{
    my $self = shift;
    my $sta_ok;
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
        my $submit = $config{'RESUBMIT'};
        my $resubmit_queued = 1;
        if (defined $config{'RESUBMIT_QUEUED'}) {
            $resubmit_queued = $config{'RESUBMIT_QUEUED'};
        }

        if ($resubmit_queued < 0) {
            die "error: RESUBMIT_QUEUED, if set, should be >= 0";
        }

        my $queuedJobs = 1 > $resubmit_queued ? 1 : $resubmit_queued;

        if ($queuedJobs >= 2 && $submit > 0) {
            print "warning: RESUBMIT_QUEUED is >= 2 which sets RESUBMIT to 0";
            my $newresubmit= 0;
            `./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit`;
        }

        while ($queuedJobs > 0) {

            my $jobname = "$config{'CASE'}.run";
            $self->addDependentJob($jobname);

            # do we add the short-term archiver to the dependency queue?
            if($config{'DOUT_S'} eq 'TRUE')
            {
                my $jobname = "$config{'CASE'}.st_archive";
                $self->addDependentJob($jobname);
            }

            $queuedJobs-=1;
        }
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
# Base class getSubmitArguments.  If we need to have submit arguments to qsub, 
# then this method will pull them out and set them up.  
#==============================================================================
sub getSubmitArguments()
{
    my $self = shift;
    # We need the script name and the dependent job id.
    my $scriptname = shift;
    my $dependentjobid = shift;

	# Get a BatchMaker instance, we need its instance data. 
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
    	    # If the argName isn't defined, just use the argflag, there;s
            # no data to replace.  
    	    elsif(defined $argFlag && ! defined $argName)
    	    {
    	        $submitargs .= " $argFlag";
    	    }
    	}
	}

    # If we have a dependent job id, we need to get the depend string
    # for this particular setup, and add it to the submit arguments. 
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
	# just to make sure, if we're here, we should be returning the
    # base class BatchUtils
    else    
    {
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
    my $batchtype = $batchsystems[0]->getAttribute('type');
    return $batchtype;
}


#==============================================================================
# Mira/ALCF specific BatchUtils class, since the workflow for ALCF has to be 
# completely different. 
# Current workflow: 
# Run the CESM run on Mira or Cetus.  When done, ssh over to tukeylogin1 and submit 
# the short-term archive run.  If we need to continue and resubmit, we will then 
# ssh back to either Mira or Cetus and resubmit the CESM run.  
#==============================================================================
package Batch::BatchUtils_mira;
use base qw( Batch::BatchUtils );
use Data::Dumper;
use Cwd;

#==============================================================================
# Overridden submitJobs() method for Mira. 
# For ALCF, we really only want this method to submit the CESM run. 
# The short-term archiver and resubmission will be handled elsewhere. 
#==============================================================================
sub submitJobs()
{
    my $self = shift;
    my $sta_ok = shift;
    my $depjobid = shift;
    
    my %depqueue = %{$self->{dependencyqueue}};

    # Get the first job sequence number. 
    my $firstjobseqnum = (sort {$a <=> $b } keys %depqueue)[0];
    # we get the first job array reference. 
    my $firstjobarray = $depqueue{$firstjobseqnum};
    # Get the first job name. 
    my $firstjobname = $$firstjobarray[0];

    # submit the CESM run, and nothing else. 
    $depjobid = $self->submitSingleJob($firstjobname, $depjobid, 0, $sta_ok);
}

#==============================================================================
# ALCF-specific single job submission. 
# The trick with ALCF is when we need to know which machine to ssh back to 
# to resubmit the CESM run.
# So, write a 'workflowhostfile' which contains the hostname we need to ssh back to 
# to resubmit the run. 
# mira submitSingleJob
#==============================================================================
sub submitSingleJob()
{
    my $self = shift;
    my $scriptname = shift;
    my $dependentJobId = shift;
    my $islastjob = shift;
    my $sta_ok = shift;
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
    #$self->SUPER::submitSingleJob($scriptname, $dependentJobId, $islastjob, $sta_ok);
    my %config = %{$self->{'caseconfig'}};
      
    my $dependarg = '';
    my $submitargs = '';
    $submitargs = $self->getSubmitArguments($scriptname, $dependentJobId);
    if(! defined $submitargs && length($submitargs <= 0))
    {
        $submitargs = '';
    }

    if(defined $sta_ok)
    {
        $submitargs .= " --env sta_ok=TRUE ";   
    }
    if(defined $islastjob)
    {
        $submitargs .= " --env islastjob=TRUE ";
    }
    
    print "Submitting CESM job script $scriptname\n";
    my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname";
    print "Runcmd: $runcmd\n";
    
    my $output;
    
    eval {
        open(my $RUN, "-|", $runcmd) // die " job submission failed, $!";
        $output = <$RUN>;
        close $RUN or die "job submission failed; |$?|, |$!|";
    };

    my $exitstatus = ($?>>8);
    if($exitstatus != 0)
    {
        print "job submission failed\n";
        exit(1);
    }
    chomp $output;
    return undef;
}
#==============================================================================
# Mira-specific doResubmit call.  If this is called from the cesm run, then we 
# have to ssh over to tukey and run the short-term archiver. 
# If called from the short-term archiver, 
#==============================================================================
sub doResubmit()
{
    my $self = shift;
    my $islastjob = shift;
    my $resubmit = shift;
    my $scriptname = shift;
    my $sta_ok = shift;

    my %config = %{$self->{'caseconfig'}};
    if(defined $sta_ok)
    {
        $ENV{'sta_ok'} = 'TRUE';
    }
    else
    {
        delete $ENV{'sta_ok'};
    }
    
    #If we're NOT doing short-term archiving, and we need to resubmit, then we need to resubmit JUST the run.  
    if($scriptname =~ /run/ && $config{'RESUBMIT'} > 0 && $config{'CONTINUE_RUN'} eq 'TRUE' && $config{'DOUT_S'} eq 'FALSE')
    {
        chdir $config{'CASEROOT'};
        my $submitargs = $self->getSubmitArguments($scriptname);
        #if($config{'RESUBMIT'} > 0)
        if($islastjob)
        {
            $submitargs .= " --env islastjob=TRUE";
        } 
         
        if(defined $sta_ok)
        {
            $submitargs .= " --env sta_ok=TRUE";
        }
        
        my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $scriptname";
        
        qx($runcmd) or die "coult not exec command $runcmd, $!";
        my $newresubmit = $config{'RESUBMIT'} - 1;
        #my $owd = getcwd;
        #chdir $config{'CASEROOT'};
        `./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit`;
        if($?)
        {
            print "could not execute ./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit\n";
        }
        #chdir $owd;

    }

    # If we ARE doing short-term archiving and we aren't resubmitting, then 
    # just run the short-term archiver 
    if($scriptname =~ /run/ && $config{'DOUT_S'} eq 'TRUE' && $config{'RESUBMIT'} == 0)
    {
        chdir $config{'CASEROOT'};
        my $starchivescript = $scriptname;
        $starchivescript =~ s/run/st_archive/g;
        my $submitargs = $self->getSubmitArguments($starchivescript);
        
        my $submitstuff = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $starchivescript";
        
        my $runcmd = "ssh tukeylogin1 $submitstuff";
    
        qx($runcmd) or die " could not exec cmd $runcmd, $! $?";
        
    }

    
    # If we're post run and we need to run the short-term archiver AND resubmit, then run the short-term archiver
    # on tukey
    if($scriptname =~ /run/ && $config{'RESUBMIT'} > 0 && $config{'CONTINUE_RUN'} eq 'TRUE' && $config{'DOUT_S'} eq 'TRUE')
    {
        chdir $config{'CASEROOT'};
        my $starchivescript = $scriptname;
        $starchivescript =~ s/run/st_archive/g;
        
        my $submitargs = $self->getSubmitArguments($starchivescript);
        if($config{'RESUBMIT'} > 0 && $config{'CONTINUE_RUN'} eq 'TRUE')
        {
            $submitargs .= " --env islastjob=TRUE ";
        }
        if(defined $sta_ok)
        {
            $submitargs .= " --env sta_ok=TRUE";
        }
            
        my $submitstuff = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $starchivescript";
        #my $cmd = "ssh tukeylogin1 qsub  -A CESM_Atmos -t 60 -n 1 -q default --mode script ./$starchivescript";
        my $runcmd = "ssh tukeylogin1 $submitstuff";
        qx($runcmd) or die "could not exec cmd $runcmd, $!";
    }

	# If we're being called by the short-term archiver, and we actually need to resubmit
	# something, then ssh from the tukey compute nodes to tukeylogin1, then ssh back to 
    # either mira or cetuslac1, and resubmit the run. 
    if($scriptname =~ /archive/ && $islastjob eq 'TRUE' && $resubmit > 0)
    {
        chdir $config{'CASEROOT'};
        my $newresubmit = $config{'RESUBMIT'} - 1;
        my $owd = getcwd;
        chdir $config{'CASEROOT'};
        `./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit`;
        chdir $owd;
        
        my $runscript = $scriptname;
        $runscript =~ s/st_archive/run/g;
        
        my $submitargs = $self->getSubmitArguments($runscript);
    
        if($config{'RESUBMIT'} > 0 && $config{'CONTINUE_RUN'} eq 'TRUE')
        {
            $submitargs .= " --env islastjob=TRUE";
        }
        if(defined $sta_ok)
        {
            $submitargs .= " --env sta_ok=TRUE";
        }
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
        if($?)
        {
            print "could not execute runcmd $runcmd, $! $?\n";
            exit(1);
        }
        $newresubmit = $config{'RESUBMIT'} - 1;
        #$owd = getcwd;
        #chdir $config{'CASEROOT'};
        `./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit`;
        if($?)
        {
            print "could not execute ./xmlchange -file env_run.xml -id RESUBMIT -val $newresubmit\n";
        }
        #chdir $owd;
        
    }
    
}

#==============================================================================
# If we need arguments when we submit a job, this is where we will add them.
# This is really specific to Mira/ALFC, no center requires qsub submit arguments. 
#==============================================================================
sub getSubmitArguments()
{
    my $self = shift;
    
    # We need the script name and the dependent job id. 
    my $scriptname = shift;
    my $dependentjobid = shift;
    my $sta_ok = shift;

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

	# Get the dependent job id if necessary..
    if(defined $dependentjobid)
    {
        my $dependArg = $self->getDependString($dependentjobid);
        $submitargs .= " $dependArg";
    }
    # Need to add the --cwd argument to the submit args if this is the st_archive script
	# so the archive script will run out of the CASEROOT
    if(defined $scriptname && $scriptname =~ /archive/)
    {
        $submitargs .= " --cwd $self->{'caseroot'} ";
    }
        
    return $submitargs;
}

#==============================================================================
package Batch::BatchUtils_cetus;
use base qw( Batch::BatchUtils_mira );
#==============================================================================
# Red herring method so the factory will work. 
#==============================================================================
sub _test()
{
    my $self = shift;
}

1;
