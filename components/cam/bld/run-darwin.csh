#! /bin/csh -f

#-----------------------------------------------------------------------
## Run script for running on a Macintosh OS-X platform (Darwin)
## using the Absoft IBM XLF/XLC compilers. This runs a low resolution
## T5 Eulerian Spectral case in serial mode.
#-----------------------------------------------------------------------
##

## Do our best to get sufficient stack memory
limit stacksize unlimited

## netCDF stuff
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = $HOME/cam_trunk

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA    $HOME/inputdata

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: startup, continue, or branch.
## $stop_n is the number of timesteps to integrate (units depends on stop_option value)
set dyn          = "eul"
set ocn          = "dom"
set case         = camrun
if ( $dyn != "eul" ) set case = "$case.$dyn"
if ( $ocn != "dom" ) set case = "$case.$ocn"
set runtype      = startup
set stop_n      = 10

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = ~/runs/
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld
set res          = "8x16"
if ( $dyn == "fv" ) set res = "10x15"

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## build exec
if ( ! -x $blddir/cam )then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    $cfgdir/configure -test -res $res -ocn $ocn -dyn $dyn -debug  || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    make -j4 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -test -case $case -runtype $runtype \
 -namelist "&camexp stop_option='nsteps', stop_n=$stop_n, nhtfrq=5/"  || echo "build-namelist failed" && exit 1

## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
mv $blddir/*in .
echo "running CAM in $rundir"
$blddir/cam                 || echo "CAM run failed" && exit 1

exit 0
