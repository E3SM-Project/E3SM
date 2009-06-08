#!/bin/csh -f

# this scripts builds pio, timing, and testpio
# edit the "USER SETTINGS" section
# run this script interactively

# ------ USER SETTINGS -----
set srcdir = `pwd`
set clean = true
set timing = true
# --------------------------

setenv NETCDF_PATH /soft/apps/netcdf-3.6.2
setenv PNETCDF_PATH /soft/apps/parallel-netcdf-1.0.2
setenv FC /bgsys/drivers/ppcfloor/comm/bin/mpixlf90_r
setenv CC /bgsys/drivers/ppcfloor/comm/bin/mpixlc_r

set confopt = ""
if ($timing =~ true) then
  set confopt = "$confopt --enable-timing"
  echo "setting confopt = $confopt"
endif

cd ../pio
./configure --disable-mct MPIF90="$FC" CC="$CC" ${confopt}

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

