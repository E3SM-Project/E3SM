#! /bin/tcsh -f
#
##=======================================================================
##
##  run-pc.csh
##
##  Generic batch submission script for PC-linux using PBS.  
##
##-----------------------------------------------------------------------
## Batch options for machine with PBS batch system.
## Usage for Lahey compiler (default): 
##   qsub run-pc.csh
## Usage for pgf90 compilers:
##   env BLD_PGI=true qsub run-pc.csh
## Usage for Intel compilers:
##   env BLD_INTEL=true qsub run-pc.csh
##-----------------------------------------------------------------------
##
## This is an example script to build and run the default CAM configuration
## on a linux system.  The default configuration is 1.9x2.5L26, Finite-Volume 
## dynamics, CLM2 land model, CICE ice model, and CAM data ocean model.  
## Script will request 1 node to run MPI-only, with 8 MPI tasks.
##
## Name of the queue (CHANGE THIS if needed)
#PBS -q long
## Number of nodes & procs/node - must have ecc memory (CHANGE THIS if needed)
#PBS -l nodes=1:ppn=8
## output file base name
#PBS -N run-pc
## Put standard error and standard out in same file
#PBS -j oe
## Export all Environment variables
#PBS -V
## End of options
##=======================================================================

##extract number of tasks from batch environment
set ntasks = `wc -l $PBS_NODEFILE`
set ntasks = $ntasks[1]

set OS = `uname -s`;
switch ( $OS )
  case Linux:
     if ( ! $?PBS_JOBID ) then
       echo "${0}: ERROR::  This batch script must be submitted via PBS";
       echo "${0}:          on a Linux machine\!";
       exit;
     else
       echo "${0}: Running CAM on Linux using PBS";
     endif

     setenv PGI /usr/local/pgi
     setenv LAHEY /usr/local/lf95

     setenv INTEL /usr/local/intel-cluster
     setenv LD_LIBRARY_PATH \
       ${PGI}/linux86/lib:${LAHEY}/lib64:/cluster/torque/lib:${INTEL}/cc/11.0.074/lib/intel64:${INTEL}/fc/11.0.074/lib/intel64:${LD_LIBRARY_PATH}
     setenv P4_GLOBMEMSIZE 500000000

     if ( $?BLD_INTEL ) then
       set netcdf = /usr/local/netcdf-intel
       set mpich = /usr/local/mpich-intel
       setenv PATH ${INTEL}/fc/11.0.074/bin/intel64:${INTEL}/cc/11.0.074/bin/intel64:${mpich}/bin:${PATH}
       ${INTEL}/intel-login-script.csh
       set cfg_string = "-fc ifort "
     else if ( $?BLD_PGI ) then
       set netcdf = /usr/local/netcdf-pgi
       set mpich = /usr/local/mpich-pgi
       setenv PATH ${PGI}/linux86/bin:${mpich}/bin:${PATH}
       set cfg_string = ""
     else
       set netcdf = /usr/local/netcdf-gcc-lf95
       set mpich = /usr/local/mpich-lf95
       setenv PATH ${LAHEY}/bin:${mpich}/bin:${PATH}
       set cfg_string = "-fc lf95 "
     endif
     setenv INC_NETCDF ${netcdf}/include
     setenv LIB_NETCDF ${netcdf}/lib
     setenv INC_MPI ${mpich}/include
     setenv LIB_MPI ${mpich}/lib
     breaksw;
  default:
    echo "${0}: This script meant for running CAM on Linux machines";    exit;
endsw

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /fs/cgd/...

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA     /fs/cgd/csm/inputdata

## LOGNAME - used in default settings, must be set if not available
## setenv LOGNAME <username>
if !($?LOGNAME) then
    echo "environment variable, LOGNAME required for setting of defaults - exiting"
    exit 1
endif

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: startup, continue, or branch.
## $stop_n is the number of days to integrate (units depends on stop_option)
set case         = camrun
set runtype      = startup
set stop_n       = 1

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = /scratch/cluster/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure $cfg_string -dyn fv -hgrid 1.9x2.5 -nosmp -ntasks $ntasks    || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $rundir                      || echo "cd $rundir failed" && exit 1
$cfgdir/build-namelist -s -config $blddir/config_cache.xml -case $case -runtype $runtype \
 -namelist "&camexp stop_option='ndays', stop_n=$stop_n /"  || echo "build-namelist failed" && exit 1

## Run CAM
echo "running CAM in $rundir"
mpiexec -n $ntasks $blddir/cam             || echo "CAM run failed" && exit 1

exit 0
