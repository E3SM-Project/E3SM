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
use Cwd;
use Data::Dumper;
use File::Basename;
use Exporter qw(import);
use XML::LibXML;
require Batch::BatchMaker;
use lib '.';
use Log::Log4perl qw(get_logger);
my $logger;

BEGIN{
    $logger = get_logger();
}

#==============================================================================
# Base class constructor.  required args are the case name, caseroot, cime root,
# compiler, machine, machine root directory, the mpi library,
# get the paths to the config_machines and config_batch xml files, and figure out
# the batch system type.
#==============================================================================
sub new
{
    my ($class, %params) = @_;
    my $self = {
	case		=> $params{'case'}		|| undef,
	caseconfig	=> $params{'caseconfig'}	|| undef,
	caseroot	=> $params{'caseroot'}		|| undef,
	cimeroot	=> $params{'cimeroot'}		|| undef,
	compiler	=> $params{'compiler'}		|| undef,
	machine		=> $params{'machine'}		|| undef,
	machroot	=> $params{'machroot'}		|| undef,
	mpilib		=> $params{'mpilib'}		|| undef,
	job => $params{job}
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
        $logger->logdie ("Could not determine batch system type for machine $self->{machine}");
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
	$logger->logdie ("Could not determine batch system type for machine $self->{machine}");
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
	$logger->logdie ("could not find depend string for this batch system type");
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
	    $logger->logdie ("could not find job id pattern for batch system type $self->{'batchtype'}");
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
	$logger->logdie (" could not ascertain dependent job id... aborting");
    }
    return $jobid;
}


#==============================================================================
# Submit jobs in a dependency chain.
#==============================================================================
sub submitJobs()
{
    my $self = shift;
    my $scriptname = shift;
    my $depjobid = undef;

    my %depqueue = %{$self->{dependencyqueue}};
    my $lastjobseqnum = (sort {$b <=> $a } keys %depqueue)[0];
    foreach my $jobnum(sort keys %depqueue)
    {
	foreach my $jobname(@{$depqueue{$jobnum}})
	{
            $logger->info("jobname: $jobname");
            $logger->debug( "lastjobseqnum $lastjobseqnum");
	    my $islastjob = 0;
	    $islastjob = 1 if ($jobnum == $lastjobseqnum);
	    $depjobid = $self->submitSingleJob($jobname, $depjobid, $islastjob);
	}
    }
    $logger->debug("in submitJobs");
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

    $logger->info("Submitting job script: $scriptname");
    #my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname $sta_argument";
    chdir $config{'CASEROOT'};
    my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname ";

    $logger->info(": $runcmd");
    my $output;

    eval {
	open (my $RUN, "-|", $runcmd) or $logger->logdie ("job submission failed, $!");
	$output = <$RUN>;
	close $RUN or $logger->logdie( "job submission failed: |$?|, |$!|");
    };
    my $exitstatus = ($?>>8);
    if($exitstatus != 0)
    {
	$logger->logdie("job submission failed $?");
    }

    chomp $output;

    my $jobid = $self->getJobID($output);
    $logger->debug( "Job ID: $jobid");
    return $jobid;
}

sub _decrementResubmitCounter()
{
    my ($self,$config) = @_;
    my $newresubmit;
    if(defined $config->{RESUBMIT}){
	$newresubmit = $config->{'RESUBMIT'} - 1;
    }else{
	$logger->logdie("RESUBMIT not defined in \$config");
    }
    my $owd = getcwd;
    chdir $config->{'CASEROOT'};
    if($config->{COMP_RUN_BARRIERS} ne "TRUE")
    {
	`./xmlchange -noecho CONTINUE_RUN=TRUE`;
    }else{
	$logger->warn("NOT changing CONTINUE_RUN since COMP_RUN_BARRIERS is on")
    }
    `./xmlchange -noecho RESUBMIT=$newresubmit`;
    if($?)
    {
	$logger->logdie( "could not execute ./xmlchange RESUBMIT=$newresubmit");
    }

}
#==============================================================================
# Base class doResubmit
# Check to see if the next set of jobs needs to be submitted.
#======================================================

sub doResubmit()
{
    #my ($self, $islastjob, $resubmit, $scriptname, $sta_ok) = @_;
    my ($self, $scriptname) = @_;

    $logger->info( "resubmitting jobs...");
    $self->dependencyCheck($scriptname);
    $self->submitJobs($scriptname);
    $self->_decrementResubmitCounter($self->{caseconfig});

}

#==============================================================================
# If we need to resubmit a job + post-run jobs, this subroutine will check the
# env*.xml variables to see which jobs need to be resubmitted.
# For now, we are only handling runs and the short-term archiver.
#==============================================================================
sub dependencyCheck()
{
    my $self = shift;
    my $scriptname = shift;;
    my %config = %{$self->{'caseconfig'}};

    $self->{dependencyqueue} = undef;
    # we always want to run the test or run again..
    if(-e "case.test")
    {
	my $jobname = "case.test";
	$self->addDependentJob($jobname);
    }
    else
    {
	my $jobname = "case.run";
	$self->addDependentJob($jobname);
    }

    # do we add the short-term archiver to the dependency queue?
    if($config{'DOUT_S'} eq 'TRUE')
    {
	my $jobname = "case.st_archive";
	$self->addDependentJob($jobname);
    }
	if($config{'DOUT_L_MS'} eq 'TRUE')
	{
		my $jobname = "case.lt_archive";
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
# Base class getSubmitArguments.  If we need to have submit arguments to qsub,
# then this method will pull them out and set them up.
#==============================================================================
sub getSubmitArguments()
{
    my $self = shift;
    # We need the script name and the dependent job id.
    my $scriptname = shift;
    my $dependentjobid = shift;
    $scriptname =~ /\w+\.(\w+)$/;
    $self->{job} = $1;
    $logger->debug(" scriptname: $scriptname job $self->{job}");
    # Get BatchMaker instance, we need its instance data.
    my $batchmaker = Batch::BatchFactory::getBatchMaker( caseroot => $self->{caseroot},
							 cimeroot => $self->{cimeroot},
							 case => $self->{case},
							 mpilib => $self->{mpilib},
							 machroot => $self->{machroot},
							 machine => $self->{machine},
							 compiler => $self->{compiler},
                                                          job => $self->{job} );


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
		    # some machines dont use the project flag
    	            $logger->warn("$argName not defined! ...");
    	        }
    	        else
    	        {
		    # if argFlag ends in an = sign dont put any space before field
		    if($argFlag =~ /=$/){
			$submitargs .= " $argFlag$field";
		    }else{
			$submitargs .= " $argFlag $field";
		    }
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

sub getBatchUtils
{
    my (%params) = @_;

    # We need a machine to be defined
    my $machine = $params{'machine'};
    if(!defined $machine)
    {
	$logger->logdie ("BatchUtilsFactory: machine must be defined!");
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
        $logger->logdie ("Could not determine batch system type for machine $machine");
    }
    my $batchtype = $batchsystems[0]->getAttribute('type');
    return $batchtype;
}


#==============================================================================
# Mira/ALCF specific BatchUtils class, since the workflow for ALCF has to be
# completely different.
# Current workflow:
# Run on Mira or Cetus.  When done, ssh over to cooleylogin1 and submit
# the short-term archive run.  If we need to continue and resubmit, we will then
# ssh back to either Mira or Cetus and resubmit the run.
#==============================================================================
package Batch::BatchUtils_mira;
use base qw( Batch::BatchUtils );

use Cwd;

#==============================================================================
# Overridden submitJobs() method for Mira.
# For ALCF, we really only want this method to submit the run.
# The short-term archiver and resubmission will be handled elsewhere.
#==============================================================================
sub submitJobs()
{
    my ($self,$scriptname) = @_;

    my %depqueue = %{$self->{dependencyqueue}};

    # Get the first job sequence number.
    my $firstjobseqnum = (sort {$a <=> $b } keys %depqueue)[0];
    # we get the first job array reference.
    my $firstjobarray = $depqueue{$firstjobseqnum};
    # Get the first job name.
    my $firstjobname = $$firstjobarray[0];

    # submit the run, and nothing else.
    $self->submitSingleJob($firstjobname);
}

#==============================================================================
# ALCF-specific single job submission.
# The trick with ALCF is when we need to know which machine to ssh back to
# to resubmit the run.
# So, write a 'workflowhostfile' which contains the hostname we need to ssh back to
# to resubmit the run.
# mira submitSingleJob
#==============================================================================
sub submitSingleJob()
{
    my ($self, $scriptname) = @_;

    my $workflowhostfile = "./workflowhostfile";
    if(! -e $workflowhostfile)
    {
        open (my $W, ">", $workflowhostfile) or $logger->logdie ("could not open workflow host file, $!");
        my $host = (defined $ENV{HOST})? $ENV{HOST}: $ENV{HOSTNAME};
        print $W "$host\n";
        close $W;
        $logger->info("Setting workflow host $host");
    }
    my %config = %{$self->{'caseconfig'}};

    my $dependarg = '';
    my $submitargs = '';
    $submitargs = $self->getSubmitArguments($scriptname);

    if(! defined $submitargs && length($submitargs <= 0))
    {
        $submitargs = '';
    }
    $logger->info( "Submitting job script $scriptname");
    my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname ";

    $logger->info( "Submitting job $runcmd");
    $logger->debug("Runcmd: $runcmd");

    my $output;

    eval {
        open(my $RUN, "-|", $runcmd) // $logger->logdie (" job submission failed, $!");
        $output = <$RUN>;
        close $RUN or $logger->logdie ("job submission failed; |$?|, |$!|");
    };

    my $exitstatus = ($?>>8);
    if($exitstatus != 0)
    {
        $logger->logdie( "job submission failed");
        exit(1);
    }
    chomp $output;
    return undef;
}
#==============================================================================
# Mira-specific doResubmit call.  If this is called from the run, then we
# have to ssh over to cooley and run the short-term archiver.
# If called from the short-term archiver,
#==============================================================================
sub doResubmit()
{
    my ($self, $scriptname) = @_;

    my %config = %{$self->{'caseconfig'}};

    #If we're NOT doing short-term archiving, and we need to resubmit, then we need to resubmit JUST the run

    my $issta = ($scriptname =~ /archive/);

    if(! $issta && $config{'RESUBMIT'} > 0  && $config{'DOUT_S'} eq 'FALSE')
    {
        chdir $config{'CASEROOT'};
        my $submitargs = $self->getSubmitArguments($scriptname);

        my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $scriptname ";
        $logger->info("1: $runcmd");
        qx($runcmd);
        if($? != 0)
        {
            $logger->logdie( "could not execute runcmd $runcmd, $! $?");
            exit(1);
        }

        #chdir $owd;
	$self->_decrementResubmitCounter(\%config);

    }

    # If we ARE doing short-term archiving and we aren't resubmitting, then
    # just run the short-term archiver
    if(! $issta && $config{'DOUT_S'} eq 'TRUE' && $config{'RESUBMIT'} == 0)
    {
        chdir $config{'CASEROOT'};

	# replace .run or .test with .st_archive

	my($basename, $path) = fileparse($scriptname);
        my $starchivescript = "$path/$basename.st_archive";
        my $submitargs = $self->getSubmitArguments($starchivescript);

        my $submitstuff = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $starchivescript";

        my $runcmd = "ssh cooleylogin1 $submitstuff";

        $logger->info("2: $runcmd");
        qx($runcmd);
        if($? != 0)
        {
            $logger->logdie( "could not execute runcmd $runcmd, $! $?");
            exit(1);
        }

    }


    # If we're post run and we need to run the short-term archiver AND resubmit, then run the short-term archiver
    # on cooley
    if(! $issta && $config{'RESUBMIT'} > 0 && $config{'DOUT_S'} eq 'TRUE')
    {
        chdir $config{'CASEROOT'};
        my $starchivescript = $scriptname;
	if($scriptname =~ /\.run$/){
	    $starchivescript =~ s/run/st_archive/g;
	}else{
	    $starchivescript =~ s/test/st_archive/g;
        }
        my $submitargs = $self->getSubmitArguments($starchivescript);

        my $submitstuff = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $starchivescript";
        my $runcmd = "ssh cooleylogin1 $submitstuff";
        $logger->info("3: $runcmd");
        qx($runcmd);
        if($? != 0)
        {
            $logger->logdie( "could not execute runcmd $runcmd, $! $?");
            exit(1);
        }

    }

    # If we're being called by the short-term archiver, and we actually need to resubmit
    # something, then ssh from the cooley compute nodes to cooleylogin1, then ssh back to
    # either mira or cetuslac1, and resubmit the run.
    if($issta && $config{RESUBMIT} > 0)
    {
        chdir $config{'CASEROOT'};

        my $runscript = $scriptname;
	if(defined $config{TESTCASE}){
	    $runscript =~ s/st_archive/test/g;
	}else{
	    $runscript =~ s/st_archive/run/g;
	}


        my $submitargs = $self->getSubmitArguments($runscript);

        my $submitstuff = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} $runscript";
        open (my $W, "<", "./workflowhostfile" ) or $logger->logdie( "could not open workflow host file, $!");
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
        my $runcmd = "ssh cooleylogin1 ssh $runhost $submitstuff ";
    #     my $runcmd =  "ssh $runhost $submitstuff ";
         $logger->info("4: $runcmd");
        qx($runcmd);
        if($? != 0)
        {
            $logger->logdie( "could not execute runcmd $runcmd, $! $?");
            exit(1);
        }
	$self->_decrementResubmitCounter(\%config);

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
    $scriptname =~ /\.([^\.]+)$/;
    my $job = $1;
    my $batchmaker = Batch::BatchFactory::getBatchMaker( caseroot => $self->{caseroot},
							 cimeroot => $self->{cimeroot},
							 case => $self->{case},
							 mpilib => $self->{mpilib},
							 machroot => $self->{machroot},
							 machine => $self->{machine},
							 compiler => $self->{compiler},
                                                         job => $job );

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
        $logger->debug( "flag: $argFlag");
        if(defined $argName && length($argName) > 0)
        {
            # Get the actual data field from the BatchMaker class.
            my $field = $batchmaker->getField($argName);
            if(! defined $field)
            {
                $logger->logdie ("$argName not defined! Aborting...");
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

    return $submitargs;
}

#==============================================================================
# Red herring method so the factory will work.
#==============================================================================
sub _test()
{
    my $self = shift;
}

1;
