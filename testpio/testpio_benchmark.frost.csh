#!/bin/csh -f

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


setenv FC /usr/bin/blrts_xlf90
setenv CC /usr/bin/gcc

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

foreach suite (all)
     set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf --enable-timing"
     set testlist = "b01 b02 b03 b04 b05 b06 b07 b08 b09 b10 b11"

  echo "Building executable for ${suite} suite..."
  echo "Configuration options for ${suite} suite are ${confopts}"
  echo "Test list for ${suite} suite are ${testlist}"

  cd ${tstdir}
  cd ../pio
  ./configure MPIF90="\$FC" CC="\$CC" ${confopts}
  gmake clean
  cd ../timing
  cp -f ../testpio/Makefile.timing ./Makefile.timing
  sed s/'$(FOPTS)'/'$(FOPTS) -WF,-DBGL'/ <Makefile.timing> Makefile
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

    echo "$suite :: $test :: qsub -n 16 -t -q debug testpio"
    set stat = `cqsub -n 64 -t 00:50:00 -q debug -o ${fout} ./testpio < testpio_in`
    cqwait $stat
##    mpirun.lsf ./testpio >>& ${fout}

   cp ${fout} ${testpiodir}/
    set pass = `grep "completed successfully" ${fout} | wc -l`
    if ($pass > 0) then
       set tstat = "PASS"
    else
       set tstat = "FAIL"
    endif
        
    echo "${tstat} ${testname} ${suite} ${test}" >> ${outfil}

  end
end

#---------------------------

###EOF

###echo "qsub -n 16 -t -q -A $project ${testpiodir}/${testname}.sub"
###qsub -n 16 -t 50 -q prod-devel -A $project  ${testpiodir}/${testname}.sub


