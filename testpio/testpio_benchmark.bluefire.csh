#!/bin/csh -f

# this sets up a suite of preset tests on bluefire
# edit the "USER SETTINGS" section
# run this script interactively on bluefire

# ------- USER SETTINGS ----
set testname = "testpiob"
set testpiodir = `pwd`
set piodir = ${testpiodir}/..
set wrkdir = "/ptmp/$USER/${testname}"
# ---------------------------

set LID = "`date +%y%m%d-%H%M%S`"
set srcdir = ${wrkdir}/src
set tstdir = ${srcdir}/testpio
set outfil = ${testpiodir}/${testname}.out.$LID

cat >! ${testpiodir}/${testname}.sub << EOF
#!/bin/csh -f

#BSUB -n 16
#BSUB -R "span[ptile=64]"
#BSUB -q premium
#BSUB -N
#BSUB -x
#BSUB -a poe
#BSUB -o poe.stdout.%J
#BSUB -e poe.stderr.%J
#BSUB -J testpio_suite
#BSUB -W 0:50
#BSUB -P 43310039

setenv NETCDF_PATH /usr/local
setenv PNETCDF_PATH /contrib/pnetcdf

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
#     set testlist = "b01 b02 b03 b04 b05 b06 b07 b08"
     set testlist = "bb03"
#  echo "confopts \${confopts}"
#  echo "testlist \${testlist}"

  cd ${tstdir}
  cd ../pio
  ./configure \${confopts}
  gmake clean
  cd ../timing
  cp -f ../testpio/Makefile.timing ./Makefile
  gmake clean
  cd ../testpio
  gmake clean
  cd ../timing
  gmake
  cd ../pio
  gmake
  cd ../testpio
  gmake

  foreach test (\${testlist})

    set casedir = ${wrkdir}/\${suite}.\${test}

    if (! -d \${casedir}) mkdir -p \${casedir}
    cd \${casedir}

    rm -f ./testpio
    cp -f ${tstdir}/testpio ./testpio
    rm -f ./testpio_in
    cp -f ${tstdir}/testpio_in.\${test} ./testpio_in
    if (! -d none) mkdir none
    rm -r -f none/*

    set fout = ${testname}.\${suite}.\${test}.out.$LID
    rm -f \${fout}
    touch \${fout}

    mpirun.lsf ./testpio >>& \${fout}

#   cp \${fout} ${testpiodir}/
    set pass = \`grep "completed successfully" \${fout} | wc -l\`
    if (\$pass > 0) then
       set tstat = "PASS"
    else
       set tstat = "FAIL"
    endif
        
    echo "\${tstat} ${testname} \${suite} \${test}" >> ${outfil}

  end
end

#---------------------------

EOF

echo "bsub  < ${testpiodir}/${testname}.sub"
bsub  < ${testpiodir}/${testname}.sub

