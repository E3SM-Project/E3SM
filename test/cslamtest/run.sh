#!/bin/tcsh -f



set HOMME = $PWD/../..
set bdir = $HOMME/build/cslam
set exe = $HOMME/build/cslam/cslam
set src = $HOMME/test/cslamtest
set wdir = ~/scratch1/cslam
mkdir $wdir
mkdir $wdir/movies
mkdir $wdir/results


set make = 1
if ( $#argv>0 ) then
if ( $1  == build  ) then
   # for climate:
#   setenv FCFLAGS "-O2 -w90"
   setenv FCFLAGS "-O0 -g  -w90"
   cd $bdir
    ./configure NC=4 NP=5  PLEV=4 NTRAC=5 --enable-pio --enable-blas --enable-lapack \
--with-pnetcdf=$PNETCDF_PATH --with-netcdf=$NETCDF_PATH
   make -j4 depends
   make clean
   make -j4
   exit
endif
endif
if ( $make ) then
  cd $bdir
  make -j4    
  if ( $status ) exit
endif


cd $wdir
set NCPU = 16
mpirun -np $NCPU $exe < $src/input.nl


