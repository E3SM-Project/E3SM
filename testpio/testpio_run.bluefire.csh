#!/bin/csh -f

# this sets up a series of tests across pe counts on bluefire
# edit the "USER SETTINGS" section
# edit the "NAMELIST INPUT" section
# run this script interactively on bluefire

# ------ USER SETTINGS -----
set run_name = 'testpio_s1'
set testpiodir = `pwd`
set new_input_file = FALSE
#set input_file = "${run_name}_in"
set input_file = "testpio_in.pb04"
set wrkdir = "/ptmp/$USER"
set pemin = 16
set pemax = 16
set ncntmax = 20
set queue = premium
set acct = 93300014

# ---- NAMELIST INPUT ------

if (${new_input_file} == TRUE) then

cat >! ${testpiodir}/${input_file} << EOF
&io_nml
  casename    = '${run_name}:pnx:box:stride=1'
  nx_global   = 360
  ny_global   = 240
  iofmt       = 'pnc'
  rearr       = 'none'
  nprocsIO    = -1
  stride      = 1
  base        = 0
  maxiter     = 10
  dir         = './none/'
  num_aggregator = 1
  DebugLevel  = 0
  compdof_input = 'namelist'
  compdof_output = 'none'
/
&compdof_nml
  nblksppe = 1
  grdorder = 'xyz'
  grddecomp = 'xy'
  gdx = 0
  gdy = 0
  gdz = 0
  blkorder = 'xyz'
  blkdecomp1 = 'xy'
  blkdecomp2 = ''
  bdx = 0
  bdy = 0
  bdz = 0
/
&prof_inparm
  profile_disable = .false.
  profile_barrier = .false.
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

#BSUB -n $pecnt
###BSUB -R "span[ptile=32]"
#BSUB -q $queue
#BSUB -N
#BSUB -x
#BSUB -a poe
#BSUB -o poe.stdout.%J
#BSUB -e poe.stderr.%J
#BSUB -J testpio
#BSUB -W 0:10
#BSUB -P $acct

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

mpirun.lsf ./testpio >& \$fout

cp \$fout \$srcdir/

EOF
#------

   @ ncnt = $ncnt + 1
   @ pecnt = $pecnt * $mfact
#  @ pecnt = $pecnt + 1

   echo "bsub < $scrname"
   bsub < $scrname

end
