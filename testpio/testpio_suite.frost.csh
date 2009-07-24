#!/bin/csh -fx

# this sets up a suite of preset tests on intrepid
# edit the "USER SETTINGS" section
# run this script interactively on intrepid

# ------- USER SETTINGS ----
set testname = "testpios"
set testpiodir = `pwd`
set piodir = ${testpiodir}/..
set wrkdir = "/ptmp/$USER/${testname}"
set project = " "
# ---------------------------

set LID = "`date +%y%m%d-%H%M%S`"
set srcdir = ${wrkdir}/src
set tstdir = ${srcdir}/testpio
set outfil = ${testpiodir}/${testname}.out.$LID

###cat >! ${testpiodir}/${testname}.sub << EOF
###!/bin/csh -f

setenv NETCDF_PATH /contrib/bgl/netcdf-3.6.2
setenv PNETCDF_PATH /contrib/bgl/parallel-netcdf-bld121807
setenv MPI_INC -I/bgl/BlueLight/ppcfloor/bglsys/include
setenv MPI_LIB '-L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts'


setenv FC blrts_xlf90
setenv F90 blrts_xlf90
setenv CC blrts_xlc
setenv FFLAGS '-g -qfullpath -qextname=flush -qflttrap=enable:underflow'

# --------------------------

if (! -d ${srcdir}) mkdir -p ${srcdir}
cp -r -f $piodir/pio $srcdir/
cp -r -f $piodir/mct $srcdir/
cp -r -f $piodir/timing $srcdir/
cp -r -f $piodir/testpio $srcdir/

if (-e ${wrkdir}/wr01.dof.txt) then
  rm -f ${wrkdir}/wr01.dof.txt
endif

touch ${outfil}

foreach suite (snet pnet mpiio all ant)
  if (${suite} =~ "snet") then
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08"
  else if (${suite} =~ "pnet") then
     set testlist = "pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08"
  else if (${suite} =~ "mpiio") then
     set testlist = "bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08"
  else if (${suite} =~ "all") then
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08 pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08 bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08 wr01 rd01"
  else if (${suite} =~ "ant") then
     set testlist = "sn02 sb02 pn02 pb02 bn02 bb02"
  else
     echo "suite ${suite} not supported"
     exit -2
  endif

  $tstdir/testpio_benchmark.frost.csh $suite $testlist
  
end

#---------------------------

###EOF

###echo "qsub -n 16 -t -q -A $project ${testpiodir}/${testname}.sub"
###qsub -n 16 -t 50 -q prod-devel -A $project  ${testpiodir}/${testname}.sub


