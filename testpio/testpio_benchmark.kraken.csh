#!/bin/csh -f

# this sets up a suite of preset tests on bluefire
# edit the "USER SETTINGS" section
# run this script interactively on bluefire

# ------- USER SETTINGS ----
set testname = "testpiob"
set testpiodir = `pwd`
set piodir = ${testpiodir}/..
set wrkdir = "/lustre/scratch/$USER/${testname}"
# ---------------------------

set LID = "`date +%y%m%d-%H%M%S`"
set srcdir = ${wrkdir}/src
set tstdir = ${srcdir}/testpio
set outfil = ${testpiodir}/${testname}.out.$LID

cat >! ${testpiodir}/${testname}.sub << EOF
#!/bin/csh -f

###PBS -A XXX
#PBS -N testpio_benchmark
#PBS -q batch
#PBS -l size=64
#PBS -l feature=1gbpercore
#PBS -l walltime=03:00:00
#PBS -j oe
#PBS -S /bin/csh -V

if (-e /opt/modules/default/init/csh) then
  source /opt/modules/default/init/csh
  module load   xtpe-quadcore
  module switch xt-mpt xt-mpt/3.1.2 # 3.0.2    is default on 2008-sep-03
  module switch pgi    pgi/7.1.6    # 7.1.6    is default on 2008-sep-03
  module load   netcdf/3.6.2        # 3.6.2    is default on 2008-sep-03
  module load   xt-craypat
  module load   p-netcdf/1.0.3
endif
setenv NETCDF_PATH \$NETCDF_DIR
setenv PNETCDF_PATH \$PNETCDF_DIR
setenv FC ftn
#setenv FFLAGS "-i4 -target=linux -gopt -Mextend -byteswapio -Mflushz -Kieee -Ktrap=fp"

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
#     set testlist = "b01 b02 b03 b04 b05 b06 b07 b08 b09 b10 b11"
     set testlist = "b09"

#  echo "confopts \${confopts}"
#  echo "testlist \${testlist}"

  cd ${tstdir}
  cd ../pio
  ./configure MPIF90="\$FC" --disable-mpi2 --enable-filesystem-hints=lustre --host=Linux \${confopts}
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

    rm testpio+pat
    pat_build -g io -o ./testpio+pat ./testpio
    aprun -n 64 ./testpio+pat >>& \${fout}

    cp \${fout} ${testpiodir}/
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

echo "qsub ${testpiodir}/${testname}.sub"
qsub ${testpiodir}/${testname}.sub

