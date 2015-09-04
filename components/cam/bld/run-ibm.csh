#! /usr/bin/csh -f

##=======================================================================
##
##  run-ibm.csh
##
##  Generic batch submission script for IBM-AIX using LSF.  
##
##-----------------------------------------------------------------------
## Usage for xlf compiler (default): 
##   bsub < run-ibm.csh
##-----------------------------------------------------------------------
##
## This is an example script to build and run the default CAM configuration
## on an IBM SP.  The default configuration is 1.9x2.5L26, Finite-Volume 
## dynamics, CLM2 land model, CICE ice model, and CAM data ocean model.  
## Script will request 1 node to run with 16 MPI tasks and 4 SMP threads/task.
##
## Setting LSF options for batch queue submission.
#BSUB -a poe                     # use LSF openmp elim
#BSUB -x                         # exclusive use of node (not_shared)
#BSUB -n 16                      # bluefire setting
#BSUB -R "span[ptile=16]"        # bluefire setting
#BSUB -o out.%J                  # output filename
#BSUB -e out.%J                  # error filename
#BSUB -q regular                 # queue
#BSUB -W 0:10                    # wall clock limit
# #BSUB -P <please_set>            # account number 

##extract number of tasks from batch environment
set ntasks = `echo $LSB_HOSTS | wc -w`

# should be set equal to (CPUs-per-node / tasks_per_node)
setenv OMP_NUM_THREADS 4

## suggestions from Jim Edwards 07/08
setenv XLSMPOPTS "stack=256000000"
setenv OMP_DYNAMIC false
setenv AIXTHREAD_SCOPE S
setenv MALLOCMULTIHEAP true
setenv MP_USE_BULK_XFER yes
setenv MP_LABELIO yes

## Do our best to get sufficient stack memory
limit stacksize unlimited

## netCDF stuff
source /contrib/Modules/3.2.6/init/csh
module load netcdf/4.1.3_seq

setenv INC_NETCDF $(NETCDF)/include
setenv LIB_NETCDF $(NETCDF)/lib

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /fis/cgd/...

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA    /fs/cgd/csm/inputdata

## LOGNAME - used in default settings, must be set if not available
## setenv LOGNAME <username>
if !($?LOGNAME) then
    echo "LOGNAME not available for setting of defaults - setting must be added to this script"
    exit 1
endif

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: startup, continue, or branch.
## $stop_n is the number of days to integrate (units depend on stop_option)
set case         = camrun
set runtype      = startup
set stop_n       = 1

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = /ptmp/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure -dyn fv -hgrid 1.9x2.5 -spmd -smp -ntasks $ntasks -nthreads $OMP_NUM_THREADS   || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $rundir                      || echo "cd $rundir failed" && exit 1
$cfgdir/build-namelist -s -config $blddir/config_cache.xml -case $case -runtype $runtype \
 -namelist "&camexp stop_option='ndays', stop_n=$stop_n /"  || echo "build-namelist failed" && exit 1

## Run CAM - use 'mpirun.lsf' on bluefire
echo "running CAM in $rundir"

#bluefire
mpirun.lsf /usr/local/bin/hybrid_launch $blddir/cam    || echo "CAM run failed" && exit 1

exit 0
