#!/bin/csh -f

# this sets up a suite of preset tests on bluefire
# edit the "USER SETTINGS" section
# run this script interactively on bluefire

# ------- USER SETTINGS ----
set testname = "testpios"
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
#PBS -N testpio_suite
#PBS -q debug
#PBS -l size=16
#PBS -l walltime=00:50:00
#PBS -j oe
#PBS -S /bin/csh -V

if (-e /opt/modules/default/init/csh) then
  source /opt/modules/default/init/csh
 module load   xtpe-quadcore
# module swap   xt-mpt xt-mpt/3.0.2 # 3.0.2    is default on 2008-sep-03
  module switch pgi    pgi/7.1.6    # 7.1.6    is default on 2008-sep-03
  module load   netcdf/3.6.2        # 3.6.2    is default on 2008-sep-03
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

foreach suite (snet pnet mpiio all ant)
  if (\${suite} =~ "snet") then
     set confopts = "--disable-mct --disable-pnetcdf --disable-mpiio --enable-netcdf --enable-timing"
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08"
  else if (\${suite} =~ "pnet") then
     set confopts = "--disable-mct --enable-pnetcdf --disable-mpiio --disable-netcdf --enable-timing"
     set testlist = "pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08"
  else if (\${suite} =~ "mpiio") then
     set confopts = "--disable-mct --disable-pnetcdf --enable-mpiio --disable-netcdf --enable-timing"
     set testlist = "bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08"
  else if (\${suite} =~ "all") then
     set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf --enable-timing"
     set testlist = "sn01 sn02 sn03 sb01 sb02 sb03 sb04 sb05 sb06 sb07 sb08 pn01 pn02 pn03 pb01 pb02 pb03 pb04 pb05 pb06 pb07 pb08 bn01 bn02 bn03 bb01 bb02 bb03 bb04 bb05 bb06 bb07 bb08 wr01 rd01"
  else if (\${suite} =~ "ant") then
     set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf --disable-timing"
     set testlist = "sn02 sb02 pn02 pb02 bn02 bb02"
  else
     echo "suite \${suite} not supported"
     exit -2
  endif

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

    aprun -n 16 ./testpio >>& \${fout}

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

echo "qsub ${testpiodir}/${testname}.sub"
qsub ${testpiodir}/${testname}.sub

