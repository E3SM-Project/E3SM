#!/bin/csh -f

# this scripts builds pio, timing, and testpio
# edit the "USER SETTINGS" section
# run this script interactively

# ------ USER SETTINGS -----
set srcdir = `pwd`
set clean = true
set timing = true
# --------------------------


if (-e /opt/modules/default/init/csh) then
  source /opt/modules/default/init/csh
  module load   xtpe-quadcore
# module swap   xt-mpt xt-mpt/3.0.2 # 3.0.2    is default on 2008-sep-03
  module switch pgi    pgi/7.1.6    # 7.1.6    is default on 2008-sep-03
  module load   netcdf/3.6.2        # 3.6.2    is default on 2008-sep-03
  module load   pnetcdf/1.0.2
endif

setenv NETCDF_PATH $NETCDF_DIR
setenv PNETCDF_PATH $PNETCDF_DIR
setenv FC ftn
#setenv FFLAGS "-i4 -target=linux -gopt -Mextend -byteswapio -Mflushz -Kieee -Ktrap=fp"

set confopt = ""
if ($timing =~ true) then
  set confopt = "$confopt --enable-timing"
  echo "setting confopt = $confopt"
endif

cd ../pio
./configure --disable-mct MPIF90="$FC" --disable-mpi2 --enable-filesystem-hints=lustre --host=Linux $confopt
#./configure --disable-mct MPIF90="$FC" FFLAGS="$FFLAGS" --disable-mpi2 --enable-filesystem-hints=lustre --host=Linux $confopt

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

