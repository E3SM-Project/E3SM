#!/usr/bin/env perl 
package Task::TaskMaker;

use strict;
use warnings;
use Data::Dumper;
use Cwd;
use POSIX qw(ceil floor);
use Exporter qw(import);

#==============================================================================
# Class constructor.  The only required parameter is the caseroot.  
# We have an array of layout strings in a specified order.  
# We then get the configuration for the case, then pull the values from the configuration 
# object using the layout strings array. 
#==============================================================================
sub new
{
	my ($class, %params) = @_;
	if(! defined $params{'caseroot'})
	{
		die "TaskMaker.pm: the caseroot must be set as an argument!";
	}
	
	my $self = {
		caseroot => $params{'caseroot'} || undef,
		removedeadtasks => $params{'removedeadtasks'} || undef,
	};
	bless $self, $class;

	# Create a set of strings needed to pull layout information out of the config 
	# object. 
    my @layoutstrings = qw/ COMP_CPL NTASKS_CPL NTHRDS_CPL ROOTPE_CPL PSTRID_CPL
                     COMP_ATM NTASKS_ATM NTHRDS_ATM ROOTPE_ATM PSTRID_ATM NINST_ATM
                     COMP_LND NTASKS_LND NTHRDS_LND ROOTPE_LND PSTRID_LND NINST_LND
                     COMP_ROF NTASKS_ROF NTHRDS_ROF ROOTPE_ROF PSTRID_ROF NINST_ROF
                     COMP_ICE NTASKS_ICE NTHRDS_ICE ROOTPE_ICE PSTRID_ICE NINST_ICE
                     COMP_OCN NTASKS_OCN NTHRDS_OCN ROOTPE_OCN PSTRID_OCN NINST_OCN
                     COMP_GLC NTASKS_GLC NTHRDS_GLC ROOTPE_GLC PSTRID_GLC NINST_GLC
                     COMP_WAV NTASKS_WAV NTHRDS_WAV ROOTPE_WAV PSTRID_WAV NINST_WAV
					 MAX_TASKS_PER_NODE PES_PER_NODE PIO_NUMTASKS PIO_ASYNC_INTERFACE /;
	$self->{'layoutstrings'} = \@layoutstrings;
	# Either the config was passed in, otherwise pull it from within the caseroot
	my %config;
	if( defined $params{'config'})
	{
		%config = $params{'config'};
		$self->{'config'} = $params{'config'};
	}
	else
	{
 	    my $toolsdir = $self->{'caseroot'} . "/Tools";
	    require ConfigCase;
    	my $cfgref = ConfigCase->new("$self->{'caseroot'}/Tools/config_definition.xml", "env_build.xml");
    	%config = $cfgref->getAllResolved();
		$self->{'config'} = \%config;
	}
	
	# Use the layout strings to pull the layout information out of the config
	foreach my $layoutstring(@{$self->{'layoutstrings'}})
	{
		$self->{$layoutstring} = $config{$layoutstring};
	}
	
  $self->{'EXEROOT'} = $config{'EXEROOT'};
  $self->{'DEFAULT_RUN_EXE_TEMPLATE_STR'} = '__DEFAULT_RUN_EXE__';
	
	# set the compiler, and MAX_TASKS_PER_NODE
	$self->{'COMPILER'} = $config{'COMPILER'};
	$self->{'MAX_TASKS_PER_NODE'} = 1 if $self->{'MAX_TASKS_PER_NODE'} < 1;
	
	# Set up arrays with the comps, tasks, threads, root PE, # instances, and pstrids 
	my @mcomps= ( $config{'COMP_CPL'}, $config{'COMP_ATM'}, $config{'COMP_LND'}, $config{'COMP_ROF'}, $config{'COMP_ICE'}, $config{'COMP_OCN'}, 
				       $config{'COMP_GLC'}, $config{'COMP_WAV'} );
	$self->{'mcomps'} = \@mcomps;
	
	my @ntasks = ( $config{'NTASKS_CPL'}, $config{'NTASKS_ATM'}, $config{'NTASKS_LND'}, $config{'NTASKS_ROF'}, $config{'NTASKS_ICE'}, $config{'NTASKS_OCN'}, 
				       $self->{'NTASKS_GLC'}, $self->{'NTASKS_WAV'} );
	$self->{'ntasks'} = \@ntasks;

	my @nthrds = ( $config{'NTHRDS_CPL'}, $config{'NTHRDS_ATM'}, $config{'NTHRDS_LND'}, $config{'NTHRDS_ROF'}, $config{'NTHRDS_ICE'}, $config{'NTHRDS_OCN'}, 
				       $config{'NTHRDS_GLC'}, $config{'NTHRDS_WAV'} );
	$self->{'nthrds'} = \@nthrds;

	my @rootpe = ( $config{'ROOTPE_CPL'}, $config{'ROOTPE_ATM'}, $config{'ROOTPE_LND'}, $config{'ROOTPE_ROF'}, $config{'ROOTPE_ICE'}, $config{'ROOTPE_OCN'}, 
				       $config{'ROOTPE_GLC'}, $config{'ROOTPE_WAV'} );
	$self->{'rootpe'} = \@rootpe;

	my @ninst = ( 1, $config{'NINST_ATM'}, $config{'NINST_LND'}, $config{'NINST_ROF'}, $config{'NINST_ICE'}, $config{'NINST_OCN'}, 
				       $config{'NINST_GLC'}, $config{'NINST_WAV'} );
	$self->{'ninst'} = \@ninst;
	my @pstrid = ( $config{'PSTRID_CPL'}, $config{'PSTRID_ATM'}, $config{'PSTRID_LND'}, $config{'PSTRID_ROF'}, $config{'PSTRID_ICE'}, $config{'PSTRID_OCN'}, 
				       $config{'PSTRID_GLC'}, $config{'PSTRID_WAV'} );
	$self->{'pstrid'} = \@pstrid;

	# At this point, we have only one method to compute the necessary values.  
	$self->_computeValues();
	return $self;

}

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

#==============================================================================
# get a new ConfigCase object, resolve the case values. 
#==============================================================================
sub _getConfig
{
	my $self = shift;

	my $toolsdir = $self->{'caseroot'} . "/Tools";
	require ConfigCase;
	
	my $config = ConfigCase->new("$self->{'caseroot'}/Tools/config_definition.xml", "$self->{'caseroot'}/env_build.xml");

	$config->getAllResolved();
	return $config;
}

#==============================================================================
# Compute the necessary numbers..
#==============================================================================
sub _computeValues
{

	my $self = shift;
	my $totaltasks = 0;
	my @ntasks = @{$self->{ntasks}};
	my @nthrds = @{$self->{nthrds}};
	my @rootpe = @{$self->{rootpe}};
	my @pstrid = @{$self->{pstrid}};

	for(my $i = 0; $i <= $#ntasks; $i++)
	{
		my $n = $ntasks[$i]; 
		my $t = $nthrds[$i]; 
		my $r = $rootpe[$i]; 
		my $p = $pstrid[$i]; 
		my $tt = $r + ($n - 1) * $p + 1;
		$totaltasks = $tt if $tt > $totaltasks;
	}
	
	# Check if we need to add pio's tasks to the total task count 
	if($self->{PIO_ASYNC_INTERFACE} eq "TRUE")
	{
		if($self->{PIO_NUMTASKS} > 0)
		{
			$totaltasks += $self->{PIO_NUMTASKS};
		}
		else
		{
			$totaltasks += $self->{PES_PER_NODE};
		}
	}
	
	# Compute max threads for each mpi task
	my @maxt;
	# first initialize maxt, max threads for each task
	for(my $i = 0; $i < $totaltasks; $i++)
	{
		$maxt[$i] = 0;
	}
	
	# calculate the max threads for each component
	for(my $i = 0; $i < $#ntasks; $i++)
	{
		my $n = $ntasks[$i];
		my $t = $nthrds[$i];
		my $r = $rootpe[$i];
		my $p = $pstrid[$i];
		
		my $c2 = 0;
        # Should be $c2 < $n
		#while($c2 > $n)
		while($c2 < $n)
		{
			my $s = $r + $c2 * $p;
			if($t > $maxt[$s]) { $maxt[$s] = $t; }
			$c2 = $c2 + 1;
		}
	}
	# remove tasks with zero threads if requested
	if(defined $self->{removedeadtasks})
	{
		my $alltasks = $totaltasks;
		for(my $c1 = 0; $c1 < $alltasks; $c1++)
		{
			if($c1 < $totaltasks && $maxt[$c1] < 1)
			{
				for( my $c3 = $c1; $c3 < $totaltasks - 1; $c3++)
				{
					$maxt[$c3] = $maxt[$c3+1];
				}
				$maxt[$totaltasks] = 0;
				$totaltasks = $totaltasks - 1;
			}
		}
	}

	# compute min/max threads over all mpi tasks and sum threads
	# reset maxt values from zero to one after checking for min values 
	# but before checking for max and summing..
	my $minthreads = $maxt[0];
	my $maxthreads = $maxt[0];

	my @sumthreads;
	$sumthreads[0] = 0;
	for(my $c1 = 1; $c1 < $totaltasks; $c1++)
	{
		if($maxt[$c1] < $minthreads) { $minthreads = $maxt[$c1];} 	
		if($maxt[$c1] < 1)           { $maxt[$c1] = 1;} 	
		if($maxt[$c1] > $maxthreads) { $maxthreads = $maxt[$c1] ;} 	
		$sumthreads[$c1] = $sumthreads[($c1-1)] + $maxt[($c1-1)];
	}
	$self->{'sumthreads'} = \@sumthreads;
	
    # Compute task and thread settings for batch commands 
	my $fullsum = 0;
	my $sum = $maxt[0];
	my $taskgeom = "(0";
	my $threadgeom = " $maxt[0]";
	my $taskcount = 1;
	my $maxtaskcount = $totaltasks;
	my $threadcount = $maxt[0];
	my $maxthreadcount = $maxt[0];
	my $aprun = "";
	my $pbsrs = "";

	my ($taskpernode, $nodecnt);
  my $totalnodecount = 0;
  my $maxtotalnodecount = 0;
	for (my $c1=1; $c1 < $totaltasks; $c1++)
	{
		$sum = $sum + $maxt[$c1];
	
		if($maxt[$c1] > $self->{'MAX_TASKS_PER_NODE'})
		{
			$fullsum += $self->{'MAX_TASKS_PER_NODE'};
			$sum = $maxt[$c1];
			$taskgeom .= $taskgeom.")($c1";
		}
		else
		{
			$taskgeom = $taskgeom.",$c1";
		}

		$threadgeom = $threadgeom . ":$maxt[$c1]";
		if($maxt[$c1]  != $threadcount)
		{
			#print "self Max_TASKS_PER_NODE: $self->{'MAX_TASKS_PER_NODE'}\n";
			#print "threadcount $threadcount\n";
			if(defined $self->{'PES_PER_NODE'}){
			    $taskpernode = min($self->{'PES_PER_NODE'},$self->{'MAX_TASKS_PER_NODE'} / $threadcount);
			}else{
			    $taskpernode = $self->{'MAX_TASKS_PER_NODE'} / $threadcount;
			}
			$taskpernode = ($taskpernode > $taskcount) ? $taskcount : $taskpernode;
			$aprun = $aprun . " -n $taskcount  -N $taskpernode -d $threadcount $self->{'DEFAULT_RUN_EXE_TEMPLATE_STR'} :";
			my $nodecount = $taskcount / $taskpernode;
      $totalnodecount += ceil($nodecount);
			$pbsrs = $pbsrs . "${nodecount}:ncpus=$self->{'MAX_TASKS_PER_NODE'}:mpiprocs=${taskpernode}:ompthreads=${threadcount}:model=";
			$threadcount = $maxt[$c1];
			$maxthreadcount = max($maxthreadcount, $maxt[$c1]);
			$taskcount = 1;
    }
		else
		{
			$taskcount += 1;
		}
	}

  my $maxtaskpernode = 1;
  $maxtaskpernode = $self->{'MAX_TASKS_PER_NODE'}/$maxthreadcount;
  if(defined $self->{'PES_PER_NODE'}){
    $maxtaskpernode = min($self->{'PES_PER_NODE'},$maxtaskpernode);
  }
  $maxtaskpernode = ($maxtaskpernode > $maxtaskcount) ? $maxtaskcount : $maxtaskpernode;
  $maxtotalnodecount = ceil($maxtaskcount/$maxtaskpernode);
	
	$fullsum = $fullsum + $sum;
	$self->{'fullsum'} = $fullsum;
	$taskgeom = $taskgeom.")";
	$self->{'taskgeom'} = $taskgeom;
	$taskpernode = min($self->{'PES_PER_NODE'},$self->{'MAX_TASKS_PER_NODE'} / $threadcount);
	$taskpernode = ($taskpernode > $taskcount) ? $taskcount : $taskpernode;
  $totalnodecount += ceil($taskcount/$taskpernode);
	if($self->{'COMPILER'} eq "intel" && $taskpernode > 1)
	{
		my $taskpernuma = ceil($taskpernode / 2);
		$aprun .= " -S $taskpernuma -cc numa_node ";
	}
	$aprun  .= " -n $taskcount -N $taskpernode -d $threadcount $self->{'DEFAULT_RUN_EXE_TEMPLATE_STR'} ";

	# add all the calculated numbers as instance data. 
	$self->{'totaltasks'} = $totaltasks;
	$self->{'taskpernode'} = $maxtaskpernode;
    $self->{'taskpernuma'}  = ceil($taskpernode / 2);
	$self->{'maxthreads'} = $maxthreads;
	$self->{'minthreads'} = $minthreads;
	$self->{'taskgeom'} = $taskgeom;
	$self->{'threadgeom'} = $threadgeom;
	$self->{'taskcount'} = $maxtaskcount;
	$self->{'threadcount'} = $maxthreadcount;
	$self->{'aprun'}  .= $aprun;
	$self->{'optnodecount'} = $totalnodecount;
	$self->{'nodecount'} = $maxtotalnodecount;

	$self->{'pbsrs'} = $pbsrs . "$self->{nodecount}:ncpus=$self->{'MAX_TASKS_PER_NODE'}:mpiprocs=${taskpernode}:ompthreads=${threadcount}:model=";

	# calculate ptile..  
	my $ptile = $self->{'MAX_TASKS_PER_NODE'} / 2;
	if($self->{'maxthreads'} > 1)
	{
		$ptile = floor($self->{'MAX_TASKS_PER_NODE'} / $self->{'maxthreads'});
	}
	$self->{'ptile'} = $ptile;
}

#=============================================================================
# has optimized mpirun command?
#=============================================================================
sub hasOptMpirun
{
  my($self, $mpirunexe) = @_;
  return (defined $self->{$mpirunexe});
}

#=============================================================================
# get optimized mpirun command
#=============================================================================
sub optMpirun
{
  my ($self, $mpirunexe, $cesmexe) = @_;
  my $optruncmd = $self->{$mpirunexe};
  $optruncmd =~ s/$self->{'DEFAULT_RUN_EXE_TEMPLATE_STR'}/$cesmexe/g;
  return $optruncmd;
}

#==============================================================================
# get the full PE sum
#==============================================================================
sub sumPES
{
	my $self = shift;
	return $self->{'fullsum'};
}

#==============================================================================
# get the sum only 
#==============================================================================
sub sumOnly
{
	my $self = shift;
	return $self->{'fullsum'};
}

#==============================================================================
# get the total tasks
#==============================================================================
sub sumTasks
{
	my $self = shift;
	return $self->{'totaltasks'};
}

#==============================================================================
# get the max # of threads
#==============================================================================
sub maxThreads
{
	my $self = shift;
	return $self->{'maxthreads'};
}

#==============================================================================
# get the task geometry string
#==============================================================================
sub taskGeometry
{
	my $self = shift;
	return $self->{'taskgeom'};
}

#==============================================================================
# get the thread geometry string
#==============================================================================
sub threadGeometry()
{
	my $self = shift;
	return $self->{'threadgeom'};
}

#==============================================================================
# get the task count.
#==============================================================================
sub taskCount()
{
	my $self = shift;
	return $self->{'taskcount'};
}
#==============================================================================
# get the thread geometry string
#==============================================================================
sub threadCount()
{
	my $self = shift;
	return $self->{'threadcount'};
}
#==============================================================================
# get the node count
#==============================================================================
sub nodeCount()
{
	my $self = shift;
	return $self->{'nodecount'};
}

#==============================================================================
# get the optimized node count
#==============================================================================
sub optNodeCount()
{
	my $self = shift;
	return $self->{'optnodecount'};
}
#==============================================================================
# TODO remove the aprun stuff from this module
# It is being constructed properly via XML configuration files. 
#==============================================================================
sub APRun()
{
	my $self = shift;
	return $self->{'aprun'};
}

#==============================================================================
# get the PBSRS string.  This should probably be refactored a better method ???
# TODO remove...
#==============================================================================
sub PBSRS()
{
	my $self = shift;
	return $self->{'pbsrs'};
}

#==============================================================================
# get the ptile.. 
#==============================================================================
sub ptile()
{
	my $self = shift;
	return $self->{'ptile'};
}

#==============================================================================
# get the tasks per node
#==============================================================================
sub taskPerNode()
{
	my $self = shift;
	return $self->{'taskpernode'};
}

#==============================================================================
# get the tasks per numa
#==============================================================================
sub taskPerNuma()
{	
    my $self = shift;
	return $self->{'taskpernuma'};
}

#==============================================================================
# get the max tasks per node
#==============================================================================
sub maxTasksPerNode()
{
	my $self = shift;
	return $self->{'MAX_TASKS_PER_NODE'};
}

#==============================================================================
# get the pe layout document for insertion into the run script
#==============================================================================
sub document()
{
	my $self = shift;
	my $doc = '';
	$doc .=  "# ---------------------------------------- \n";
	$doc .=  "# PE Layout: \n";
	$doc .=  "#   Total number of tasks: $self->{'totaltasks'}\n";
	$doc .=  "#   Maximum threads per task: $self->{'maxthreads'}\n";
	for(my $c1 = 0; $c1 < $#{ $self->{'ntasks'}}; $c1++ )
	{
		my $n = ${$self->{'ntasks'}}[$c1];
		my $t = ${$self->{'nthrds'}}[$c1];
		my $r = ${$self->{'rootpe'}}[$c1];
		my $i = ${$self->{'ninst'}}[$c1];
		my $p = ${$self->{'pstrid'}}[$c1];
		my $tt = $r + ($n -1) * $p;
		$doc .=  "#     " . ${$self->{'mcomps'}}[$c1] . " ntasks=$n nthreads=$t rootpe=$r ninst=$i \n";
	}
		$doc .=  "#    \n";
		$doc .=  "#    total number of hw pes = $self->{'fullsum'} \n";
	for (my $c1 = 0; $c1 < $#{$self->{'ntasks'}}; $c1++)
	{	
		my $n  = ${$self->{'ntasks'}}[$c1];
		my $t  = ${$self->{'nthrds'}}[$c1];
		my $r  = ${$self->{'rootpe'}}[$c1];
		my $p  = ${$self->{'pstrid'}}[$c1];
		my $tt = $r + ($n -1) * $p;
		my $tm = ${$self->{'sumthreads'}}[$tt] + $t + 1;
		$doc .=  "#     " . ${$self->{'mcomps'}}[$c1] . " hw pe range ~ from " . ${$self->{'sumthreads'}}[$r] . " to $tm\n";
	}

	if($self->{'minthreads'} < 1)
	{
		$doc .=  "#   \n";
    	$doc .=  "#   WARNING there appear to be some IDLE hw pes \n";
    	$doc .=  "#   Please consider reviewing your env_mach_pes.xml file \n";
	}
	$doc .=  "# ---------------------------------------- \n";
	
	return $doc;
}
    

1;
