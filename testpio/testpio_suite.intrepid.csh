#!/bin/csh -f

# this sets up a suite of preset tests on intrepid
# edit the "USER SETTINGS" section
# run this script interactively on intrepid

# ------- USER SETTINGS ----
set testname = "testpios"
set testpiodir = `pwd`
set piodir = ${testpiodir}/..
set wrkdir = "/gpfs1/$USER/${testname}"
set project = CCESDev
# ---------------------------

set LID = "`date +%y%m%d-%H%M%S`"
set srcdir = ${wrkdir}/src
set tstdir = ${srcdir}/testpio
set outfil = ${testpiodir}/${testname}.out.$LID

###cat >! ${testpiodir}/${testname}.sub << EOF
###!/bin/csh -f

setenv NETCDF_PATH /soft/apps/netcdf-3.6.2
setenv PNETCDF_PATH /soft/apps/parallel-netcdf-1.0.2

setenv F90 /bgsys/drivers/ppcfloor/comm/bin/mpixlf90_r
setenv FC /bgsys/drivers/ppcfloor/comm/bin/mpixlf90_r
setenv CC /bgsys/drivers/ppcfloor/comm/bin/mpixlc_r
setenv FFLAGS '-qextname=flush'

#setenv MPIF90 /bgsys/drivers/ppcfloor/comm/bin/mpixlf90_r
#setenv MPICC /bgsys/drivers/ppcfloor/comm/bin/mpixlc_r

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
     set confopts = "--disable-mct --disable-pnetcdf --disable-mpiio --enable-netcdf --enable-timing"
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08"
  else if (${suite} =~ "pnet") then
     set confopts = "--disable-mct --enable-pnetcdf --disable-mpiio --disable-netcdf --enable-timing"
     set testlist = "pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08"
  else if (${suite} =~ "mpiio") then
     set confopts = "--disable-mct --disable-pnetcdf --enable-mpiio --disable-netcdf --enable-timing"
     set testlist = "bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08"
  else if (${suite} =~ "all") then
     set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf --enable-timing"
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08 pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08 bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08 wr01 rd01"
  else if (${suite} =~ "ant") then
     set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf --disable-timing"
     set testlist = "sn02 sb02 pn02 pb02 bn02 bb02"
  else
     echo "suite ${suite} not supported"
     exit -2
  endif

  echo "Building executable for ${suite} suite..."
  echo "Configuration options for ${suite} suite are ${confopts}"
  echo "Test list for ${suite} suite are ${testlist}"

  cd ${tstdir}
  cd ../pio
  ./configure MPIF90="$FC" CC="$CC" --build powerpc-bgp-linux --host powerpc64-suse-linux ${confopts}
  gmake clean
  cd ../timing
  cp -f ../testpio/Makefile.timing ./Makefile.timing
  sed s/'$(FOPTS)'/'$(FOPTS) -WF,-DBGP'/ <Makefile.timing> Makefile
  rm Makefile.timing
  gmake clean
  cd ../testpio
  gmake clean
  cd ../timing
  gmake
  cd ../pio
  gmake
  cd ../testpio
  gmake

  foreach test (${testlist})

    set casedir = ${wrkdir}/${suite}.${test}

    if (! -d ${casedir}) mkdir -p ${casedir}
    cd ${casedir}

    rm -f ./testpio
    cp -f ${tstdir}/testpio ./testpio
    rm -f ./testpio_in
    cp -f ${tstdir}/testpio_in.${test} ./testpio_in
    if (! -d none) mkdir none
    rm -r -f none/*

    set fout = ${testname}.${suite}.${test}.out.$LID
    rm -f ${fout}
    touch ${fout}

    echo "$suite :: $test :: qsub -n 16 -t -q -A $project testpio"
    qsub -n 16 -t 50 -q prod-devel -A $project testpio -o ${fout}
##    mpirun.lsf ./testpio >>& ${fout}

#   cp ${fout} ${testpiodir}/
#    set pass = \`grep "completed successfully" ${fout} | wc -l\`
#    if ($pass > 0) then
#       set tstat = "PASS"
#    else
#       set tstat = "FAIL"
#    endif
#        
#    echo "${tstat} ${testname} ${suite} ${test}" >> ${outfil}

  end
end

#---------------------------

###EOF

###echo "qsub -n 16 -t -q -A $project ${testpiodir}/${testname}.sub"
###qsub -n 16 -t 50 -q prod-devel -A $project  ${testpiodir}/${testname}.sub


