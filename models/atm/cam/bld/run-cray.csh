#! /bin/tcsh -f
#
##=======================================================================
##
##  run-cray.csh
##
##  Generic batch submission script for CRAY XT using PBS.  
##
##-----------------------------------------------------------------------
## Usage for PGI compilers (default): 
##   qsub run-cray.csh
##-----------------------------------------------------------------------
##
## This is an example script to build and run the default CAM configuration
## on a CRAY XT.  The default configuration is 1.9x2.5L26, Finite-Volume 
## dynamics, CLM2 land model, CICE ice model, and CAM data ocean model.  
## Script will request 8 nodes to run MPI-only, with 32 MPI tasks.
##
# Name of the queue (CHANGE THIS if needed)
#PBS -q debug
# Number of procs and walltime (CHANGE THIS if needed)
#PBS -l walltime=0:58:00,size=32
# output file base name
#PBS -N run-cray
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
#PBS -A <please_set>
# End of options
#=======================================================================

source $MODULESHOME/init/csh
module load netcdf/3.6.2

setenv INC_NETCDF ${NETCDF_DIR}/include
setenv LIB_NETCDF ${NETCDF_DIR}/lib

setenv MPICH_MAX_SHORT_MSG_SIZE 1024

## set this equal to #nodes X #ppn
set procs = 32

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /tmp/work/$LOGNAME/...

## ROOT OF CAM DATA DISTRIBUTION
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA     /lustre/widow1/proj/ccsm/inputdata

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
set wrkdir       = /lustre/widow1/scratch/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure -dyn fv -hgrid 1.9x2.5 -spmd -nosmp -ntasks $procs -fc ftn    || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $rundir                      || echo "cd $rundir failed" && exit 1
$cfgdir/build-namelist -s -config $blddir/config_cache.xml -case $case -runtype $runtype \
 -namelist "&camexp stop_option='ndays' stop_n=$stop_n  /"  || echo "build-namelist failed" && exit 1

## Run CAM
echo "running CAM in $rundir"
aprun -n $procs $blddir/cam             || echo "CAM run failed" && exit 1

exit 0
