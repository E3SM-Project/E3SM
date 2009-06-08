#!/bin/csh -f

# this sets up a series of tests across pe counts on bluefire
# edit the "USER SETTINGS" section
# edit the "NAMELIST INPUT" section
# run this script interactively on bluefire

# ------ USER SETTINGS -----
set run_name = 'testpio_s1'
set testpiodir = `pwd`
set new_input_file = TRUE
set input_file = "${run_name}_in"
#set input_file = "testpio_in.pb02"
set wrkdir = "/lustre/scratch/$USER"
set pemin = 800
set pemax = 800
set ncntmax = 20
set queue = batch
set acct = "TG-ATM090011"

# ---- NAMELIST INPUT ------

if (${new_input_file} == TRUE) then

cat >! ${testpiodir}/${input_file} << EOF
&io_nml
  casename    = '${run_name}:pnx:box:stride=1'
  nx_global   = 3600
  ny_global   = 2400
  nz_global   = 40
  iofmt       = 'bin'
  rearr       = 'box'
  nprocsIO    = -1
  stride      = 10
  base        = 0
  maxiter     = 10
  dir         = './none/'
  num_aggregator = 8
  DebugLevel  = 0
  compdof_input = 'namelist'
  compdof_output = 'none'
  iodof_input = 'namelist'
/
&compdof_nml
  nblksppe = 1
  grdorder = 'xyz'
  grddecomp = 'setblk'
  gdx = 180
  gdy = 60
  gdz = 40
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
  gdy = 60
  gdz = 10
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

# -----------------------

set pecnt = $pemin
set ncnt = 1
set mfact = 2
while (($pecnt <= $pemax) && ($ncnt <= $ncntmax))
   set scrname = ${run_name}.${pecnt}.sub

#------
cat >! ${scrname} << EOF
#!/bin/csh -f

#PBS -A ${acct}
#PBS -N testpio
#PBS -q ${queue}
#PBS -l size=${pecnt}
#PBS -l walltime=00:20:00
#PBS -j oe
#PBS -S /bin/csh -V

set srcdir = ${testpiodir}
set casedir = "${wrkdir}/${run_name}.${pecnt}"
set LID = "`date +%y%m%d-%H%M%S`"
set fout = "${run_name}.out.${pecnt}.\$LID"

if (! -d \$casedir) mkdir \$casedir
cd \$casedir
rm -f ./testpio
cp -f \$srcdir/testpio ./testpio
rm -f ./testpio_in
cp -f \$srcdir/${input_file} ./testpio_in
if (! -d none) mkdir none
rm -r -f none/*

aprun -n ${pecnt} ./testpio >& \$fout

cp \$fout \$srcdir/

EOF
#------

   @ ncnt = $ncnt + 1
   @ pecnt = $pecnt * $mfact
#  @ pecnt = $pecnt + 1

   echo "qsub $scrname"
   qsub $scrname

end
