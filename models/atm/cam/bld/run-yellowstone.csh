#! /bin/csh -f

##=======================================================================
##
##  run-yellowston.csh
##
##  Generic batch submission script for yellowstone using LSF.  
##
##-----------------------------------------------------------------------
## Usage for intel compiler (default): 
##   ./run_yellowstone.csh build
##   bsub < run-yellowstone.csh
##-----------------------------------------------------------------------
##
## This is an example script to build and run the default CAM configuration
## on Yellowstone.  The default configuration is 1.9x2.5L26, Finite-Volume 
## dynamics, CLM2 land model, CICE ice model, and CAM data ocean model.  
## Script will request 1 node to run with 16 MPI tasks and 2 SMP threads/task.
##
## Setting LSF options for batch queue submission.
#BSUB -a poe                     # use LSF openmp elim
#BSUB -x                         # exclusive use of node (not_shared)
#BSUB -n 16                      # yellowstone setting
#BSUB -R "span[ptile=16]"        # yellowstone setting
#BSUB -o out.%J                  # output filename
#BSUB -e out.%J                  # error filename
#BSUB -q small                   # queue
#BSUB -W 0:10                    # wall clock limit
# #BSUB -P <please_set>            # account number 

if ($#argv == 1 && $1 != "build") then
  echo "Invalid command line option: $1 "
  exit
endif

##set number of tasks
set ntasks = 16

# should be set equal to (CPUs-per-node / tasks_per_node)
setenv OMP_NUM_THREADS 2

## Do our best to get sufficient stack memory
limit stacksize unlimited

## Recommended OMP setting from Jim Edwards
setenv OMP_STACKSIZE 256M

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /fis/cgd/...

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA    $CESMDATAROOT/inputdata/

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
set wrkdir       = /glade/scratch/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam && $1 == "build") then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure -dyn fv -hgrid 1.9x2.5 -spmd -smp -ntasks $ntasks -nthreads $OMP_NUM_THREADS -cc mpicc -fc mpif90 -fc_type intel  || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
# Test to see if build was done
else if (! -e $blddir/cam) then
  echo "ERROR: $blddir/cam does not exists.  You need to run run-yellowstone.sh build."
  exit
endif

## Create the namelist
cd $rundir                      || echo "cd $rundir failed" && exit 1
$cfgdir/build-namelist -s -config $blddir/config_cache.xml -case $case -runtype $runtype \
 -namelist "&camexp stop_option='ndays', stop_n=$stop_n /"  || echo "build-namelist failed" && exit 1

if ( $1 != "build") then
  ## Run CAM - use 'mpirun.lsf' on yellowstone
  echo "running CAM in $rundir"
  
  #yellowstone
  mpirun.lsf $blddir/cam    || echo "CAM run failed" && exit 1
endif

exit 0
