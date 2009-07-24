#!/bin/csh -f

# this sets up a suite of preset tests on intrepid
# edit the "USER SETTINGS" section
# run this script interactively on intrepid

# ------- USER SETTINGS ----
set run_name = "testpior"
set testpiodir = `pwd`
set new_input_file = TRUE
set piodir = ${testpiodir}/..
set wrkdir = "/ptmp/$USER/${run_name}"
set input_file = "testpio_in"
set project = " "
set nodes =  64
set mode = vn
set queue = debug
# ---------------------------

set LID = "`date +%y%m%d-%H%M%S`"
set srcdir = ${wrkdir}/src
set tstdir = ${srcdir}/testpio
set outfil = ${testpiodir}/${run_name}.out.$LID

###cat >! ${testpiodir}/${run_name}.sub << EOF
###!/bin/csh -f

setenv NETCDF_PATH /contrib/bgl/netcdf-3.6.2
setenv PNETCDF_PATH /contrib/bgl/parallel-netcdf-bld121807
setenv MPI_INC -I/bgl/BlueLight/ppcfloor/bglsys/include
setenv MPI_LIB '-L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts'


setenv FC /usr/bin/blrts_xlf90
setenv CC /usr/bin/gcc

# ---- NAMELIST INPUT ------

echo "hello #1"
if (${new_input_file} == TRUE) then

cat >! ${testpiodir}/${input_file} << EOF
&io_nml
  casename    = '${run_name}:pnx:box:stride=1'
  nx_global   = 3600
  ny_global   = 2400
  nz_global   = 1
  iofmt       = 'pnc'
  rearr       = 'box'
  nprocsIO    = -1
  stride      = 8
  base        = 0
  maxiter     = 10
  dir         = './none/'
  num_aggregator = 4
  DebugLevel  = 0
  compdof_input = 'namelist'
  compdof_output = 'none'
  iodof_input = 'namelist'
/
&compdof_nml
  nblksppe = 1
  grdorder = 'xyz'
  grddecomp = 'setblk'
  gdx = 450
  gdy = 150
  gdz = 1
  blkorder = 'xyz'
  blkdecomp1 = 'xy'
  blkdecomp2 = ''
  bdx = 0
  bdy = 0
  bdz = 0
/
&iodof_nml
  nblksppe = 1
  grdorder = 'xyz'
  grddecomp = 'setblk'
  gdx = 3600
  gdy = 150
  gdz = 1
  blkorder = 'xyz'
  blkdecomp1 = 'yz'
  blkdecomp2 = ''
  bdx = 0
  bdy = 0
  bdz = 0
/
&prof_inparm
  profile_disable = .false.
  profile_barrier = .true.
  profile_single_file = .false.
  profile_depth_limit = 10
  profile_detail_limit = 0
/
EOF

else
#   echo "${testpiodir} ${input_file}"
   if (! -e ${testpiodir}/${input_file}) then
      echo "testpio input file does not exist ${testpiodir}/${input_file}"
      exit -9
   endif
endif
echo "hello #2"

# -----------------------

#   set scrname = ${run_name}.${nodes}.sub
   set scrname = ${run_name}.sub
echo "hello #3"

#------
##cat >! ${scrname} << EOF

echo "hello #3.1"
# --------------------------

if (! -d ${srcdir}) mkdir -p ${srcdir}
cp -r -f $piodir/pio $srcdir/
cp -r -f $piodir/mct $srcdir/
cp -r -f $piodir/timing $srcdir/
cp -r -f $piodir/testpio $srcdir/

if (-e ${wrkdir}/wr01.dof.txt) then
  rm -f ${wrkdir}/wr01.dof.txt
endif

echo "hello #4"
touch ${outfil}

 set confopts = "--disable-mct --enable-pnetcdf --enable-mpiio --enable-netcdf --enable-timing"

echo "point #5"

#  cd ${testdir}
#  cd ../pio
#  ./configure MPIF90="\$FC" CC="\$CC"
#  gmake clean
#  cd ../timing
#  cp -f ../testpio/Makefile.timing ./Makefile.timing
#  sed s/'$(FOPTS)'/'$(FOPTS) -WF,-DBGL'/ <Makefile.timing> Makefile
#  rm Makefile.timing
#  gmake clean
#  cd ../testpio
#  gmake clean
#  cd ../timing
#  gmake
#  cd ../pio
#  gmake
#  cd ../testpio
#  gmake


    set casedir = ${wrkdir}/all

    if (! -d ${wrkdir}) mkdir -p ${wrkdir}
    cd ${wrkdir}
    rm -f ./testpio
    cp -f ${tstdir}/testpio ./testpio
    if (! -d none) mkdir none
    rm -r -f none/*

    rm -rf ${input_file}
    cp ${tstdir}/${input_file} ${wrkdir}/
    echo "point #8"

    set fout = ${run_name}.out.$LID
    rm -f ${fout}
    touch ${fout}
    echo "point #9"

    set stat = `cqsub -n ${nodes} -m ${mode} -t 00:50:00 -q debug -o ${fout} ./testpio < ${input_file}`
    cqwait $stat
##    mpirun.lsf ./testpio >>& ${fout}

    echo "point #10"
   cp ${fout} ${testpiodir}/
    set pass = `grep "completed successfully" ${fout} | wc -l`
    if ($pass > 0) then
       set tstat = "PASS"
    else
       set tstat = "FAIL"
    endif
        

#---------------------------

