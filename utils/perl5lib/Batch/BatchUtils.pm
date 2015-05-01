#!/usr/bin/env perl;

use strict;
use warnings;
#use diagnostics;
package Batch::BatchUtils;
use Data::Dumper;
use Cwd;
use Exporter qw(import);
use XML::LibXML;
require Batch::BatchMaker;
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
    print "BatchUtils constructor\n";
    #print Dumper $self;

    my $perl5libdir = undef;
	my $configbatch = $self->{'machroot'} . "/config_batch.xml";
	$self->{'configbatch'} = $configbatch;
    print "self configbatch $self->{'configbatch'}\n";
	my $configmachines = $self->{'machroot'} . "/config_machines.xml";
	$self->{'configmachines'} = $configmachines;
    my $casetoolsdir = $self->{'caseroot'} . "/Tools";
    push(@INC, $casetoolsdir);
    my $xml = XML::LibXML->new(no_blanks => 1);
    my $machineconfig = $xml->parse_file($configmachines);
    my $root = $machineconfig->getDocumentElement();
	
    #my @batchtypes = $root->findnodes("/batch/machine[\@MACH=\'$self->{machine}\']/batchtype");
    my @batchtypes = $root->findnodes("/config_machines/machine[\@MACH=\'$self->{machine}\']/batch_system");
    if(! @batchtypes)
    {
        die "Could not determine batch system type for machine $self->{machine}";
    }
    #$self->{'batchtype'} = $batchtypes[0]->textContent();
    $self->{'batchtype'} = $batchtypes[0]->getAttribute('name');

	$self->{dependencyqueue} = undef;
    #print Dumper $self;
	return $self;
}
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
		die "could not find depend argument for this batch system type";
	}
	my $deparg = $dependargs[0]->textContent();
    $deparg =~ s/jobid/$jobid/g;
	return $deparg;
}

# If we need arguments when we submit a job, this is where we will add them. 
sub getSubmitArguments()
{
    my $self = shift;
    my $dependentjobid = shift;
    print "in getSubmitArguments()\n";
    my $xml = XML::LibXML->new(no_blanks => 1);
    my $batchconfig = $xml->parse_file($self->{'configbatch'});
    my $root = $batchconfig->getDocumentElement();
    print "self batchtype: $self->{'batchtype'}\n";
    my @dependargs = $root->findnodes("/config_batch/batch_system[\@type=\'$self->{'batchtype'}\']/submit_args/arg");
    print "depend args: \n";
    print Dumper \@dependargs;
    my $submitargs = '';
    my $batchmaker = Batch::BatchFactory::getBatchMaker( caseroot => $self->{caseroot}, case => $self->{case}, 
                                                  mpilib => $self->{mpilib}, scriptsroot => $self->{scriptsroot},
                                                  machroot => $self->{machroot}, machine => $self->{machine},
                                                  compiler => $self->{compiler} );
    foreach my $dependarg(@dependargs)
    {
        
        my $argFlag = $dependarg->getAttribute('flag');
        my $argName = $dependarg->getAttribute('name');
        print "flag: $argFlag\n";
        if(defined $argName && length($argName) > 0)
        {
            print "name: $argName\n";
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
    print "submit args are now: $submitargs\n";
    return $submitargs;

}
# Get the job id from the output of the job submission. 
sub getJobID()
{
	my $self = shift;
	my $jobstring = shift;
	chomp $jobstring;
	#print "jobstring : |$jobstring|\n";
	my $xml = XML::LibXML->new(no_blanks => 1);
    my $batchconfig = $xml->parse_file($self->{'configbatch'});
    my $root = $batchconfig->getDocumentElement();
	#print "self batchtype:" ,  $self->{'batchtype'} , "\n";
	my @jobidpatterns = $root->findnodes("/config_batch/batch_system[\@type=\'$self->{'batchtype'}\']/jobid_pattern");
	if(!@jobidpatterns)
	{
		die "could not find job id pattern for batch system type $self->{'batchtype'}";
	}
	
	my $jobid = undef;
	my $jobidpat = $jobidpatterns[0]->textContent();
	my $pattern  = qr($jobidpat);
	#print "pattern: $pattern\n";
	if($jobstring =~ /$pattern/ )
	{
		$jobid = $1;
	}
	#print "job id is now: $jobid\n";
	return $jobid;
}


# Submit jobs in a dependency chain.
sub submitJobs()
{
	my $self = shift;
	my $depjobid = undef;
	
	my %depqueue = %{$self->{dependencyqueue}};
	my $lastjobseqnum = (sort {$b <=> $a } keys %depqueue)[0];
	#print "last dependent job number: $lastjobseqnum\n";
	#print "last dependent job: ", @{$depqueue{$lastjobseqnum}} ,  "\n";
	foreach my $jobnum(sort keys %depqueue)
	{
		foreach my $jobname(@{$depqueue{$jobnum}})
		{
			#print "job name $jobname\n";
			my $islastjob = 0;
			$islastjob = 1 if ($jobnum == $lastjobseqnum);
			$depjobid = $self->submitSingleJob($jobname, $depjobid, $islastjob);
			#print "dependent job id: $depjobid\n";
		}
	}
}

# Submit a job, given the script name. If a dependent job id is given, 
# use it as the dependency. If the islastjob flag is passed in, add the 
# islastjob=TRUE flag to the job submission.. 
sub submitSingleJob()
{
	my $self = shift;
	my $scriptname = shift;
	my $dependentJobId = shift;
	#print "in submitSingleJob, dependendJobId is: $dependentJobId\n";
	my $islastjob = shift;
	my %config = %{$self->{'caseconfig'}};
	my $dependarg = '';
	my $submitargs = '';
    $submitargs = $self->getSubmitArguments($dependentJobId);
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
	print "config batchredirect: $config{'BATCHREDIRECT'}\n";
	my $runcmd = "$config{'BATCHSUBMIT'} $submitargs $config{'BATCHREDIRECT'} ./$scriptname";
    
	print "runcmd is $runcmd\n";
	open (my $RUN, "-|", $runcmd) or die "cannot run the command \"$runcmd\"";
	#local $/ = undef;
	my $output = <$RUN>;
	chomp $output;	
	print "submit output: |$output|\n";
	
    #print "last dependent job number: $lastjobseqnum\n";
    #print "last dependent job: ", @{$depqueue{$lastjobseqnum}} ,  "\n";
	
	my $jobid = $self->getJobID($output);
	return $jobid;
}


sub doResubmit()
{
	my $self = shift;
	my %config = %{$self->{caseconfig}};
	#print "in checkforResubmit\n";
	#print "resubmitting jobs!!!\n";		
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

# If we need to resubmit a CESM job + post-run jobs, this subroutine will check the 
# env*.xml variables to see which jobs need to be resubmitted.  
sub dependencyCheck()
{
	my $self = shift;
	my %config = %{$self->{'caseconfig'}};
	# we always want to run the CESM test or run again..
	if(-e "$config{'CASE'}.test")
	{
		my $jobname = "$config{'CASE'}.test";
		$self->addDependentJob($jobname);
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
	# Do we add the tseries_generator to the dependency queue? 
	if($config{'DOUT_S_GENERATE_TSERIES'} eq 'TRUE')
	{
		my $jobname = "$config{'CASE'}.pyreshaper";
		$self->addDependentJob($jobname);
	}
	# TODO Handle submitting diagnostics here...
	
	# Do we add the long-term archiver to the dependency queue?
	if($config{'DOUT_L_MS'} eq 'TRUE')
	{
		my $jobname = "$config{'CASE'}.pyreshaper";
		$self->addDependentJob($jobname);
	}
}

# Adds a job script name or array of job script names to the dependent jobs queue. 
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

	# The jobref is a regular string scalar, make an array out of it 
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
	#print "Dependency queue now looks like: \n";
	#print Dumper $self->{dependencyqueue};
}


package BatchUtilsFactory;
use Exporter qw(import);
use XML::LibXML;
sub getBatchUtils
{
	my $machine = shift;
	my $machroot = shift;
	my $caseroot = shift;
	if(!defined $machine)
	{
		die "BatchUtilsFactory: machine must be defined!";
	}
	#$self->{'batchtype'} = getBatchSystemType($machine, $machroot, $caseroot);
	

}
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
	my $batchtype = $batchsystems[0]->textContent();
	return $batchtype;
}

package BatchUtils_lsf;
use Exporter qw(import);
sub getJobID()
{
}

sub setJobDependency()
{
}
package BatchUtils_pbs;
use Exporter qw(import);

package BatchUtils_slurm;
use Exporter qw(import);
1;
