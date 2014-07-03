#!/bin/csh
#
# testcases.csh
#
# Test the different options to the PTCLM python script.
#
# This goes through most of the different options and cases
# There are a few that are missing however specifically:
#
#     --list, --namelist, --rmold, --sitegroupname, and --debug
#
# Environment variables to set:
#
# CESM_ROOT:        To test with a separate root to CLM/CESM code set the
#                   env variable CESM_ROOT to the root directory to use.
#
# CLM_SOFF:         If set to TRUE -- stop on first failed test
#
# CLM_RETAIN_FILES: If set to FALSE -- cleanup tools build first
#
# DEBUG:            If set to TRUE -- setup cases, but do not build or run
#

set pwd=`pwd`
set mycsmdata=$HOME/inputdata
set host=`hostname`
set casedir="$pwd/myPTCLMtests.$$"
echo "Run testing for PTCLM.py on $host"

#
# Get path to root
#
if ( ! $?CESM_ROOT )then
   cd "../../../../../.."
   setenv CESM_ROOT `pwd`
   cd -
endif
setenv CCSMROOT $CESM_ROOT
if ( ! $?CLM_SOFF )then
   setenv CLM_SOFF "FALSE"
endif
if ( ! $?CLM_RETAIN_FILES )then
   setenv CLM_RETAIN_FILES "TRUE"
endif
if ( ! $?DEBUG )then
   setenv DEBUG "FALSE"
endif
#
# Machine dependent stuff
#
unset SCRATCH
if (  $host =~ ys* )then
  echo "Setup for yellowstone"
  module load netcdf/4.3.0-rc4
  module load ncl
  set parcmp=32
  set machine="yellowstone"
  set compiler="intel"
  set csmdata=/glade/p/cesm/cseg/inputdata
  set rundata="/glade/scratch/$USER"
  set netcdf=$NETCDF
  set toolsmake=""
  # Setup hostfile so mkmapdata can be run interactively
  hostname > hostfile
  setenv MP_HOSTFILE `pwd`/hostfile
  setenv MP_PROCS 1
else if ( $host =~ frankfurt* )then
  echo "Setup for frankfurt"
  set parcmp=2
  set machine="frankfurt"
  set compiler="pgi"
  set csmdata=/fs/cgd/csm/inputdata
  set rundata="/scratch/cluster/$USER"
  set netcdf=/usr/local/netcdf-pgi
  set toolsmake=""
  setenv PATH "${PATH}:/usr/bin"
else if ( $host =~ yongi* || $host =~ vpn* )then
  echo "Setup for yongi"
  set parcmp=12
  set machine="userdefined"
  set compiler="intel"
  set csmdata=/fs/cgd/csm/inputdata
  set rundata="/ptmp/$USER"
  set SCRATCH=$rundata
  set netcdf="/opt/local"
  set toolsmake="USER_FC=ifort USER_LINKER=ifort USER_CC=icc "
  set case='$CASE'
  set xmlchangeuser="./xmlchange OS=Darwin,MAX_TASKS_PER_NODE=1,MPILIB=mpi-serial,RUNDIR=/ptmp/$USER/$case/run,DIN_LOC_ROOT=$mycsmdata,COMPILER=intel,EXEROOT=/ptmp/$USER/$case,GMAKE=make"
  setenv NETCDF_PATH $netcdf
else if ( $host =~ titan* )then
  echo "Setup for titan"

  source /opt/modules/default/init/csh
  module load szip
  module load hdf5
  module load netcdf
  module load p-netcdf
  module load esmf/5.2.0rp2
  module load subversion
  module load cmake
  set netcdf=$NETCDF_PATH
  set parcmp=9
  set machine="titan"
  set compiler="pgi"
  set csmdata=/tmp/proj/ccsm/inputdata
  set rundata="/tmp/work/$USER"
  set toolsmake="USER_FC=ftn USER_CC=cc "
else
  echo "Bad host to run on: know about yellowstone, frankfurt, yong, and titan"
  exit -3
endif
alias xmlchangeuser
setenv INC_NETCDF ${netcdf}/include
setenv LIB_NETCDF ${netcdf}/lib
#
# Create or update the links to my csmdata location
#
echo "Make sure datasets are properly softlinked"
$CESM_ROOT/scripts/link_dirtree $csmdata $mycsmdata
if ( $status != 0 ) exit -1
#
# Build the tools
#
echo "Build the tools"
cd $CESM_ROOT/models/lnd/clm/tools/clm4_5/mksurfdata_map/src
if ( $CLM_RETAIN_FILES == FALSE || (! -x mksurfdata_map) && $DEBUG != TRUE )then
   gmake clean
   gmake OPT=TRUE SMP=TRUE -j $parcmp $toolsmake
   if ( $status != 0 ) exit -1
   gmake clean
endif
cd $CESM_ROOT/tools/mapping/gen_domain_files/src
if ( $CLM_RETAIN_FILES == FALSE || (! -x gen_domain) && $DEBUG != TRUE )then
   $CESM_ROOT/scripts/ccsm_utils/Machines/configure -mach $machine -compiler $compiler
   gmake clean
   gmake OPT=TRUE -j $parcmp $toolsmake
   if ( $status != 0 ) exit -1
   gmake clean
endif
cd $pwd
#
# Test the different compsets and a couple different sites
# make sure both supported compsets and flux tower sites are used
#
set caseprefix="PTCLM.$$"
set statuslog="tc.$$.status"
@ casenum = 1
echo "Write status info to $statuslog"
cat << EOF  > $statuslog
PTCLM Single-Point Simulation testing Status Log on $host


Testcase                              	          Test Status
EOF
mkdir -p $casedir
foreach mysite ( 1x1_mexicocityMEX US-UMB )
  if ( $mysite == "1x1_mexicocityMEX" || $mysite == "1x1_vancouverCAN" || $mysite == "1x1_brazil" || $mysite == "1x1_urbanc_alpha" ) then
     set suprted=TRUE
  else
     set suprted=FALSE
  endif
  if ( "$suprted" == "TRUE" ) set compsets = ( ICLM45CN ICLM45 ICLM45 )
  if ( "$suprted" != "TRUE" ) then
     if ( $mysite == "US-UMB" ) then
        set compsets = ( I_1850_CLM45 I20TRCLM45 ICLM45CN I1850CLM45BGC ICLM45BGC ICLM45BGC ICLM45BGC )
     else
        set compsets = ( ICLM45 ICLM45 ICLM45 )
     endif
  endif
  set n=0
  set finidat="none"
  foreach compset ( $compsets )
    if ( $compset == ICLM45 ) @ n = $n + 1
    set opt="--caseidprefix=$casedir/${caseprefix}_${n}"
    set opt="$opt --cesm_root $CESM_ROOT"
    if ( "$suprted" == "TRUE" ) then
      set opt="$opt --nopointdata --stdurbpt"
    else
      set opt="$opt --owritesrf --run_units=ndays --run_n=5"
    endif
    if ( $compset == I20TRCLM45 || $compset == I20TRCLM45CN ) set opt="$opt --coldstart"
    set casename="${caseprefix}_${n}_${mysite}_${compset}"
    # Use QIAN forcing on second "I" compset
    if ( $n == 2 ) then
      set opt="$opt --useQIAN --QIAN_tower_yrs"
      set casename="${casename}_QIAN"
    endif
    set case = "$casedir/$casename"
    # Use global PFT and SOIL on third (and use initial conditions from previous)
    if ( $n == 3 ) then
      if ( $finidat != "none" ) set opt="$opt --finidat $finidat"
      set opt="$opt --pftgrid --soilgrid --verbose"
    else
      set opt="$opt --quiet"
    endif
    \rm -rf $rundata/$casename
    echo "Run PTCLM for $casename options = $opt"
    set msg="$casenum $casename.PTCLM\t\t\t"
    echo    "$msg"
    echo -n "$msg" >> $statuslog
    @ casenum = $casenum + 1
    ./PTCLM.py -d $mycsmdata -m ${machine}_${compiler} -s $mysite -c $compset $opt
    if ( $status != 0 )then
       echo "FAIL $status" 		>> $statuslog
       if ( "$CLM_SOFF" == "TRUE" ) exit -1
    else
       echo "PASS"         		>> $statuslog
    endif
    cd $case
    if ( $machine == "userdefined" ) `$xmlchangeuser`
    ./xmlchange DOUT_S=FALSE
    if ( $status != 0 ) exit -1
    set msg="$casenum $casename.config\t\t\t"
    echo    "$msg"
    echo -n "$msg" >> $statuslog
    ./cesm_setup
    if ( $status != 0 )then
       echo "FAIL $status" 		>> $statuslog
       if ( "$CLM_SOFF" == "TRUE" ) exit -1
    else
       echo "PASS"         		>> $statuslog
    endif
    set msg="$casenum $casename.build\t\t\t"
    echo    "$msg"
    echo -n "$msg" >> $statuslog
    if ( $DEBUG != "TRUE" )then
       ./$casename.build
    else
       set status=1
    endif
    if ( $status != 0 )then
       echo "FAIL $status" 		>> $statuslog
       if ( "$CLM_SOFF" == "TRUE" ) exit -1
    else
       echo "PASS"         		>> $statuslog
    endif
    set msg="$casenum $casename.run\t\t\t"
    echo    "$msg"
    echo -n "$msg" >> $statuslog
    if ( $DEBUG != "TRUE" )then
       ./$casename.run
    else
       set status=1
    endif
    source ./Tools/ccsm_getenv || exit -2
    gunzip $RUNDIR/cpl.log*.gz
    set CplLogFile = `ls -1t $RUNDIR/cpl.log* | head -1`
    set basestatus = "UNDEF"
    if ( $status != 0 )then
       set basestatus = "FAIL"
    else
       set pass = `grep "SUCCESSFUL TERM" $CplLogFile | wc -l`
       if ( $pass != 0 ) then
         set basestatus = "PASS"
       else
         set basestatus = "FAIL"
       endif
    endif
    echo "$basestatus"  		>> $statuslog
    if ( "$CLM_SOFF" == "TRUE" && $basestatus == "FAIL" ) exit -1
    if ( $DEBUG != "TRUE" )then
       set finidat=`ls -1 $rundata/$casename/run/$casename.clm?.r.*.nc | tail -1`
    else
       set finidat="$rundata/$casename/run/$casename.clm2.r.0001-01-01-00000.nc"
    endif
    # Clean the build up
    #./$casename.clean_build
    if ( $compset != ICLM45 && $n > 1 ) set n = 0
    cd $pwd
  end
end
set closemsg="Successfully ran all test cases for PTCLM"
echo
echo
echo $closemsg
echo $closemsg >> $statuslog
