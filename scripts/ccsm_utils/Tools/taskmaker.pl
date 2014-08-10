#!/usr/bin/env perl

use strict;
use Cwd;

#-------------------------------------------------------------------------------
# This script derives task and thread batch settings
# based on environment variable values set in the env_configure.xml file
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# parse arg vector: select option
#-------------------------------------------------------------------------------

my $sumpeflag    = -1; # total number of pes: (mpi tasks)x(threads)
my $sumtasks     = -1; # total number of mpi tasks
my $maxthrds     = -1; # max threads over all mpi tasks
my $taskgeomflag = -1; # task geometry string for IBM
my $thrdgeomflag = -1; # thread geometry string for IBM
my $aprunflag    = -1; # aprun options for Cray XT
my $pbsrsflag    = -1; # pbs resources option for NAS (NASA pleiades)
my $nas_node_type    ; # NAS node type name (har|wes|san)
my $document     = -1; # document the layout
my $removedeadtasks = 0; # remove dead tasks or reset to 1

my ($sumpeflag, $sumtaskflag,$opt,$maxthrdflag,$pbsrsflag,$tottasks,$c1,$n,$t,$r);
if ($#ARGV  < 0 ){ # no arguments
   $taskgeomflag = 1;  # default option
}
elsif ($#ARGV eq 0 ){ # one argument
#  $opt=shift(@ARGV);
   $opt=$ARGV[0];
   if     ($opt eq "-sumonly"){
      $sumpeflag = 1;
   }
   elsif ($opt eq "-sumpes"){
      $sumpeflag = 1;
   }
   elsif ($opt eq "-sumtasks"){
      $sumtaskflag = 1;
   }
   elsif ($opt eq "-maxthrds"){
      $maxthrdflag = 1;
   }
   elsif ($opt eq "-taskgeom"){
      $taskgeomflag = 1;
   }
   elsif ($opt eq "-thrdgeom"){
      $thrdgeomflag = 1;
   }
   elsif ($opt eq "-aprun"){
      $aprunflag = 1;
   }
   elsif ($opt eq "-document"){
      $document = 1;
   }
   else {
      print "(taskmaker.pl) Usage: taskmaker.pl [-taskgeom|-thrdgeom|-sumpes|-sumtasks|-maxthrds|-aprun|-pbsrs (har|wes|san)] \n";
      exit;
   }
}
elsif ($#ARGV eq 1 ){ # two arguments
   $opt=$ARGV[0];
   if ($opt eq "-pbsrs"){
       $pbsrsflag = 1; $nas_node_type = $ARGV[1];
   }
   else {
      print "(taskmaker.pl) Usage: taskmaker.pl [-taskgeom|-thrdgeom|-sumpes|-sumtasks|-maxthrds|-aprun|-pbsrs (har|wes|san)] \n";
      exit;
   }
}
else { # more than two arguments
   print "(taskmaker.pl) Usage: taskmaker.pl [-taskgeom|-thrdgeom|-sumpes|-sumtasks|-maxthrds|-aprun|-pbsrs (har|wes|san)] \n";
   exit;
}

#-------------------------------------------------------------------------------
# parse the xml files to create xmlvars hash
#-------------------------------------------------------------------------------
my $caseroot = $ENV{CASEROOT};
my $dirs = "${caseroot}/Tools";
unshift @INC, $dirs;
require XML::Lite;

my %xmlvars = ();


my @files = <${caseroot}/*xml>;
foreach my $file (@files) {
    my $xml = XML::Lite->new( "$file" );
    my @e = $xml->elements_by_name('entry');
    while ( my $e = shift @e ) {
	my %a = $e->get_attributes();
	$xmlvars{$a{'id'}} = $a{'value'};
    }
}

#-------------------------------------------------------------------------------
# get task/thread layout data via env vars
#-------------------------------------------------------------------------------

my $COMP_CPL     = $xmlvars{'COMP_CPL'};
my $NTASKS_CPL   = $xmlvars{'NTASKS_CPL'};
my $NTHRDS_CPL   = $xmlvars{'NTHRDS_CPL'};
my $ROOTPE_CPL   = $xmlvars{'ROOTPE_CPL'};
my $PSTRID_CPL   = $xmlvars{'PSTRID_CPL'};
my $COMP_ATM     = $xmlvars{'COMP_ATM'};
my $NTASKS_ATM   = $xmlvars{'NTASKS_ATM'};
my $NTHRDS_ATM   = $xmlvars{'NTHRDS_ATM'};
my $ROOTPE_ATM   = $xmlvars{'ROOTPE_ATM'};
my $NINST_ATM    = $xmlvars{'NINST_ATM'};
my $PSTRID_ATM   = $xmlvars{'PSTRID_ATM'};
my $COMP_LND     = $xmlvars{'COMP_LND'};
my $NTASKS_LND   = $xmlvars{'NTASKS_LND'};
my $NTHRDS_LND   = $xmlvars{'NTHRDS_LND'};
my $ROOTPE_LND   = $xmlvars{'ROOTPE_LND'};
my $NINST_LND    = $xmlvars{'NINST_LND'};
my $PSTRID_LND   = $xmlvars{'PSTRID_LND'};
my $COMP_ROF     = $xmlvars{'COMP_ROF'};
my $NTASKS_ROF   = $xmlvars{'NTASKS_ROF'};
my $NTHRDS_ROF   = $xmlvars{'NTHRDS_ROF'};
my $ROOTPE_ROF   = $xmlvars{'ROOTPE_ROF'};
my $NINST_ROF    = $xmlvars{'NINST_ROF'};
my $PSTRID_ROF   = $xmlvars{'PSTRID_ROF'};
my $COMP_ICE     = $xmlvars{'COMP_ICE'};
my $NTASKS_ICE   = $xmlvars{'NTASKS_ICE'};
my $NTHRDS_ICE   = $xmlvars{'NTHRDS_ICE'};
my $ROOTPE_ICE   = $xmlvars{'ROOTPE_ICE'};
my $NINST_ICE    = $xmlvars{'NINST_ICE'};
my $PSTRID_ICE   = $xmlvars{'PSTRID_ICE'};
my $COMP_OCN     = $xmlvars{'COMP_OCN'};
my $NTASKS_OCN   = $xmlvars{'NTASKS_OCN'};
my $NTHRDS_OCN   = $xmlvars{'NTHRDS_OCN'};
my $ROOTPE_OCN   = $xmlvars{'ROOTPE_OCN'};
my $NINST_OCN    = $xmlvars{'NINST_OCN'};
my $PSTRID_OCN   = $xmlvars{'PSTRID_OCN'};
my $COMP_GLC     = $xmlvars{'COMP_GLC'};
my $NTASKS_GLC   = $xmlvars{'NTASKS_GLC'};
my $NTHRDS_GLC   = $xmlvars{'NTHRDS_GLC'};
my $ROOTPE_GLC   = $xmlvars{'ROOTPE_GLC'};
my $NINST_GLC    = $xmlvars{'NINST_GLC'};
my $PSTRID_GLC   = $xmlvars{'PSTRID_GLC'};
my $COMP_WAV     = $xmlvars{'COMP_WAV'};
my $NTASKS_WAV   = $xmlvars{'NTASKS_WAV'};
my $NTHRDS_WAV   = $xmlvars{'NTHRDS_WAV'};
my $ROOTPE_WAV   = $xmlvars{'ROOTPE_WAV'};
my $NINST_WAV    = $xmlvars{'NINST_WAV'};
my $PSTRID_WAV   = $xmlvars{'PSTRID_WAV'};
my $MAXTPN       = $xmlvars{'MAX_TASKS_PER_NODE'};
my $PESPN        = $xmlvars{'PES_PER_NODE'};
my $PIO_NUMTASKS = $xmlvars{'PIO_NUMTASKS'};
my $PIO_ASYNC_INTERFACE = $xmlvars{'PIO_ASYNC_INTERFACE'};

my $COMPILER = $xmlvars{COMPILER};

if ($MAXTPN < 1) {$MAXTPN = 1 ;}

my @mcomps = (  $COMP_CPL,   $COMP_ATM,   $COMP_LND,   $COMP_ICE,   $COMP_OCN,   $COMP_GLC,   $COMP_WAV,   $COMP_ROF);
my @ntasks = ($NTASKS_CPL, $NTASKS_ATM, $NTASKS_LND, $NTASKS_ICE, $NTASKS_OCN, $NTASKS_GLC, $NTASKS_WAV, $NTASKS_ROF);
my @nthrds = ($NTHRDS_CPL, $NTHRDS_ATM, $NTHRDS_LND, $NTHRDS_ICE, $NTHRDS_OCN, $NTHRDS_GLC, $NTHRDS_WAV, $NTHRDS_ROF);
my @rootpe = ($ROOTPE_CPL, $ROOTPE_ATM, $ROOTPE_LND, $ROOTPE_ICE, $ROOTPE_OCN, $ROOTPE_GLC, $ROOTPE_WAV, $ROOTPE_ROF);
my @ninst  = (          1,  $NINST_ATM,  $NINST_LND,  $NINST_ICE,  $NINST_OCN,  $NINST_GLC,  $NINST_WAV,  $NINST_ROF);
my @pstrid = ($PSTRID_CPL, $PSTRID_ATM, $PSTRID_LND, $PSTRID_ICE, $PSTRID_OCN, $PSTRID_GLC, $PSTRID_WAV, $PSTRID_ROF);

#print "ntasks = @ntasks \n";
#print "nthrds = @nthrds > $MAXTHREADSPERTASK\n";
#print "rootpe = @rootpe \n";
#print "ninst  = @ninst \n";
#print "pstrid = @pstrid \n";

#-------------------------------------------------------------------------------
# compute total number of mpi tasks
#-------------------------------------------------------------------------------

$tottasks = 0;
for ($c1=0; $c1 <= $#ntasks; $c1++){
    my $n = $ntasks[$c1];
    my $t = $nthrds[$c1];
    my $r = $rootpe[$c1];
    my $p = $pstrid[$c1];

    my $tt = $r + ($n - 1) * $p + 1;
    if ($tt > $tottasks) {$tottasks = $tt ;}
}
if($PIO_ASYNC_INTERFACE eq "TRUE"){
    if($PIO_NUMTASKS>0) {
	$tottasks += $PIO_NUMTASKS;
    }else{
	$tottasks += $PESPN;
    }
}
#-------------------------------------------------------------------------------
# compute max threads for each mpi task
#-------------------------------------------------------------------------------
my @maxt;
# initialize maxt, max threads for each task
for ($c1=0; $c1 < $tottasks; $c1++){
    $maxt[$c1] = 0;
}

# compute maxt array (max threads for each task)
for ($c1=0; $c1 <= $#ntasks; $c1++){
    my $n = $ntasks[$c1];
    my $t = $nthrds[$c1];
    my $r = $rootpe[$c1];
    my $p = $pstrid[$c1];

    my $c2 = 0;
    while ($c2 < $n) {
       my $s = $r + $c2 * $p;
       if ($t > $maxt[$s]) {$maxt[$s] = $t;}
       $c2 = $c2 + 1;
    }
}

# remove tasks with zero threads if requested
if ($removedeadtasks > 0) {
  my $alltasks = $tottasks;
  for ($c1=0; $c1 < $alltasks; $c1++){
    if ($c1 < $tottasks && $maxt[$c1] < 1) {
      for (my $c2=$c1; $c2 < $tottasks-1; $c2++){
        $maxt[$c2] = $maxt[$c2+1];
      }
      $maxt[$tottasks] = 0;
      $tottasks = $tottasks - 1;
    }
  }
}

# compute min/max threads over all mpi tasks and sum threads
# also reset maxt values from zero to one after checking for min values
# but before checking for max and summing
my $minthrds = $maxt[0];
my $maxthrds = $maxt[0];
my @sumt;
$sumt[0] = 0;
for ($c1=1; $c1 < $tottasks; $c1++){ 
   if ($maxt[$c1] < $minthrds) {$minthrds = $maxt[$c1] ;}
   if ($maxt[$c1] < 1) {$maxt[$c1] = 1;}
   if ($maxt[$c1] > $maxthrds) {$maxthrds = $maxt[$c1] ;}
   $sumt[$c1] = $sumt[($c1-1)] + $maxt[($c1-1)];
}

#-------------------------------------------------------------------------------
# compute task & thread settings for batch commands
#-------------------------------------------------------------------------------

my $fullsum = 0;     # sum of all tasks on all nodes
my $sum = $maxt[0];  # sum of all tasks on one node
my $taskgeom = "(0";
my $thrdgeom = " $maxt[0]";
my $taskcnt = 1;
my $thrdcnt = $maxt[0];
my $aprun = "";
my $pbsrs = "";

my ($taskpernode, $nodecnt);
for ($c1=1; $c1 < $tottasks; $c1++){     # assign each task to a node
    $sum = $sum + $maxt[$c1];
    if ($sum > $MAXTPN) {
        $fullsum = $fullsum + $MAXTPN;
        $sum = $maxt[$c1];
        $taskgeom = $taskgeom.")($c1";   # this is 1st task on a new node
    }
    else {
        $taskgeom = $taskgeom.",$c1";    # append this task to current node
    }
    $thrdgeom = $thrdgeom.":$maxt[$c1]"; # number of threads assigned to this task
    if ($maxt[$c1] != $thrdcnt) {
      $taskpernode = $MAXTPN / $thrdcnt;
      $taskpernode = ($taskpernode > $taskcnt) ? $taskcnt : $taskpernode;
      $aprun = $aprun." -n $taskcnt -N $taskpernode -d $thrdcnt \${EXEROOT}/cesm.exe :";
      $nodecnt = $taskcnt / $taskpernode ;
      $pbsrs = $pbsrs."${nodecnt}:ncpus=${MAXTPN}:mpiprocs=${taskpernode}:ompthreads=${thrdcnt}:model=${nas_node_type}+";
      $thrdcnt = $maxt[$c1];
      $taskcnt = 1;
    }
    else {
      $taskcnt = $taskcnt + 1;
    }
}
$fullsum = $fullsum + $sum;
$taskgeom = $taskgeom.")";
$taskpernode = $MAXTPN / $thrdcnt;
$taskpernode = ($taskpernode > $taskcnt) ? $taskcnt : $taskpernode;
if ($COMPILER eq "intel" && $taskpernode>1){
    my $taskpernuma = $taskpernode/2;
    $aprun .= " -S $taskpernuma -cc numa_node ";
}
$aprun .= " -n $taskcnt -N $taskpernode -d $thrdcnt \${EXEROOT}/cesm.exe";


$nodecnt = $taskcnt / $taskpernode ;
$pbsrs = $pbsrs."${nodecnt}:ncpus=${MAXTPN}:mpiprocs=${taskpernode}:ompthreads=${thrdcnt}:model=${nas_node_type}";

#print "taskgeom = $taskgeom \n";

#-------------------------------------------------------------------------------
# output what was asked for
#-------------------------------------------------------------------------------

#print "test output1 = $fullsum $tottasks $maxthrds \n";
#print "test output2 = $taskgeom \n";
#print "test output3 = $thrdgeom \n";
#print "test output4 = $aprun \n";

if    ($sumpeflag    > 0) {
    print "$fullsum";
}
elsif ($sumtaskflag  > 0) {
    print " $tottasks";
}
elsif ($maxthrdflag  > 0) {
    print "$maxthrds";
}
elsif ($taskgeomflag > 0) {
    print "$taskgeom";
}
elsif ($thrdgeomflag > 0) {
    print "$thrdgeom";
}
elsif ($aprunflag > 0) {
    print "$aprun";
}
elsif ($pbsrsflag > 0) {
    print "$pbsrs";
}
elsif ($document > 0) {
    print "# ---------------------------------------- \n";
    print "# PE LAYOUT: \n";
    print "#   total number of tasks  = $tottasks \n";
    print "#   maximum threads per task = $maxthrds \n";
  for ($c1=0; $c1 <= $#ntasks; $c1++) {
    my $n = $ntasks[$c1];
    my $t = $nthrds[$c1];
    my $r = $rootpe[$c1];
    my $i = $ninst[$c1];
    my $p = $pstrid[$c1];
    my $tt = $r + ($n - 1) * $p;
    print "#   $mcomps[$c1] ntasks=$n  nthreads=$t rootpe=$r ninst=$i \n";
  }
    print "#   \n";
    print "#   total number of hw pes = $fullsum \n";
  for ($c1=0; $c1 <= $#ntasks; $c1++) {
    my $n = $ntasks[$c1];
    my $t = $nthrds[$c1];
    my $r = $rootpe[$c1];
    my $p = $pstrid[$c1];
    my $tt = $r + ($n - 1) * $p;
    my $tm = $sumt[$tt] + $t - 1;
    print "#     $mcomps[$c1] hw pe range ~ from $sumt[$r] to $tm \n";
  }
  if ($minthrds < 1) {
    print "#   \n";
    print "#   WARNING there appear to be some IDLE hw pes \n";
    print "#   Please consider reviewing your env_mach_pes.xml file \n";
    }
    print "# ---------------------------------------- \n";
}
else {
    print "(taskmaker.pl) internal output selection error";
}

exit;

#===============================================================================
