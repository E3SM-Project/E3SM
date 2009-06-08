#!/bin/csh -f

# this scripts builds pio, timing, and testpio
# edit the "USER SETTINGS" section
# run this script interactively

# ------ USER SETTINGS -----
set srcdir = `pwd`
set clean = true
set timing = true
# --------------------------

setenv NETCDF_PATH /usr/local
setenv PNETCDF_PATH /contrib/pnetcdf

set confopt = ""
if ($timing =~ true) then
  set confopt = "$confopt --enable-timing"
  echo "setting confopt = $confopt"
endif

cd ../pio
./configure --disable-mct OPT="-O2 -g" $confopt
#./configure --disable-mct --enable-mpiio OPT="-O2 -g" $confopt
#./configure --disable-mct --enable-mpiio OPT="-O2 -g -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -qspillsize=4000 -C" $confopt
#./configure --disable-mct $confopt

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

