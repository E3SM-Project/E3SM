#!/bin/csh -f

# this scripts builds pio, timing, and testpio
# edit the "USER SETTINGS" section
# run this script interactively

# ------ USER SETTINGS -----
set srcdir = `pwd`
set clean = true
set timing = false
# --------------------------

setenv NETCDF_PATH /contrib/bgl/netcdf-3.6.2
setenv PNETCDF_PATH /contrib/bgl/parallel-netcdf-bld121807
setenv MPI_INC -I/bgl/BlueLight/ppcfloor/bglsys/include
setenv MPI_LIB '-L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts'

setenv F90 /usr/bin/blrts_xlf90
setenv FC  /usr/bin/blrts_xlf90
setenv CC  /usr/bin/gcc


set confopt = ""
if ($timing =~ true) then
  set confopt = "$confopt --enable-timing"
  echo "setting confopt = $confopt"
endif

cd ../pio
./configure --enable-pnetcdf --enable-mpiio --enable-netcdf --disable-mct MPIF90="$FC" CC="$CC" ${confopt}

if ($clean =~ true) then
  echo "cleaning first"
  cd ../timing
  cp -f ../testpio/Makefile.timing ./Makefile
  gmake clean
  cd ../pio
  gmake clean
  cd ../testpio
  gmake clean
endif

cd ../timing
cp -f ../testpio/Makefile.timing ./Makefile
gmake

cd ../pio
gmake

cd ../testpio
gmake

