#!/bin/csh -fx

set suite = $1
set n=$#argv
set testrequest = "$argv[2-$n]"

# this sets up a suite of preset tests on frost
# edit the "USER SETTINGS" section
# run this script interactively on frost

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

  if (${suite} =~ "snet") then
     set confopts = "--disable-mct --disable-pnetcdf --disable-mpiio --enable-netcdf "
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08"
  else if (${suite} =~ "pnet") then
     set confopts = "--disable-mct --enable-pnetcdf --disable-mpiio --disable-netcdf "
     set testlist = "pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08"
  else if (${suite} =~ "mpiio") then
     set confopts = "--disable-mct --disable-pnetcdf --enable-mpiio --disable-netcdf "
     set testlist = "bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08"
  else if (${suite} =~ "all") then
     set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf "
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08 pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08 bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08 wr01 rd01"
  else if (${suite} =~ "ant") then
     set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf "
     set testlist = "sn02 sb02 pn02 pb02 bn02 bb02"
  else
     echo "suite ${suite} not supported"
     exit -2
  endif


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



  echo "Building executable for ${suite} suite..."
  echo "Configuration options for ${suite} suite are ${confopts}"
  echo "Test list for ${suite} suite are ${testlist}"

  cd ${tstdir}

  perl testpio_build.pl ${confopts}

  foreach test (${testrequest})

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


