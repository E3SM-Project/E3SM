#!/bin/tcsh -f



set HOMME = $PWD/../..
set bdir = $HOMME/build/preqx
set exe = $HOMME/build/preqx/preqx
set src = $HOMME/test/cslamtest
set wdir = ~/scratch1/cslam
mkdir $wdir
mkdir $wdir/movies
mkdir $wdir/results
echo $#argv

set make = 1
if ( $#argv > 0 ) then
if ( $1 == build ) then
   # for climate:
    setenv FCFLAGS "-O2 -w90"
#   setenv FCFLAGS "-g -w90"
   cd $bdir
#    ./configure NP=4 PLEV=26 --with-arch=Darwin 
    ./configure NP=4 PLEV=26 \
--enable-blas --enable-lapack --enable-pio \
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
#cp $src/baro.nl input.nl
cp $src/baro-ne30nc4.nl input.nl

cp $src/../vcoord/cam?-26.fbin.littleendian .

set NCPU = 32
#set NCPU = 2
mpirun -np $NCPU $exe < input.nl


