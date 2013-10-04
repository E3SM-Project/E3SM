#!/usr/bin/sh
. /opt/modules/default/init/sh 

module rm PrgEnv-intel
module rm PrgEnv-pgi
module rm PrgEnv-cray 
module rm PrgEnv-gnu
module rm PrgEnv-pathscale
module rm intel
module rm pgi
module rm cray
module rm pathscale
module rm netcdf

module load PrgEnv-pgi
module switch pgi       pgi/12.5.0
module switch cray-mpich2 cray-mpich2/5.5.2
module switch cray-libsci cray-libsci/11.1.00
module load netcdf/4.1.1.0
export NETCDF_DIR=/opt/cray/netcdf/4.1.1.0/netcdf-pgi
rm cprnc
gmake  LIB_NETCDF=$NETCDF_DIR/lib INC_NETCDF=$NETCDF_DIR/include NETCDF=$NETCDF_DIR LDFLAGS="-L$NETCDF_DIR/lib -lnetcdff -lnetcdf"