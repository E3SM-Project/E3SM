#!/bin/csh

source /usr/share/Modules/init/csh
module purge
module load intel/13.0
module load openmpi-intel/1.6

setenv NETCDFROOT /projects/ccsm/yellowstone/netcdf-4.3.2-intel-13.0-openmpi-1.6
setenv PATH $NETCDFROOT/bin:$PATH
setenv LD_LIBRARY_PATH $NETCDFROOT/lib:$LD_LIBRARY_PATH
setenv NETCDF_INCLUDES $NETCDFROOT/include
setenv NETCDF_LIBS $NETCDFROOT/lib

setenv PNETCDFROOT $NETCDFROOT

setenv USER_FC ifort

gmake  LIB_NETCDF=$NETCDF_LIBS INC_NETCDF=$NETCDF_INCLUDES NETCDF=$NETCDFROOT LDFLAGS="-L$NETCDF_LIBS -lnetcdff -lnetcdf"#! /bin/csh -f
