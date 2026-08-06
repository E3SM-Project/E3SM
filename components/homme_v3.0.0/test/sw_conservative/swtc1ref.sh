#!/bin/tcsh -f
#
#SBATCH --job-name baroclinic
#SBATCH -N 4
#SBATCH --account=FY104068
#SBATCH --time=0:10:00
#
#XPBS -l nodes=4:ppn=2
#XPBS -l walltime=0:10:00
#XPBS -N swtc1
#XPBS -j oe
#XPBS -A FY081407
#XPBX -a 1758

###############################################################
# RKSSP default benchmark (used to check nothing is broken)
###############################################################
#
# NE=11, dt=360, nu=0  limiter=4  filter_freq=0 NV=4
# set smooth = 0
# set integration = runge_kutta
# set rk_stage = 3
#
# Error for cosine bell, gaussian, cyclinder, slotted cyclinter
#k= 1 12.00 days  l1,l2,linf=  0.3251433E+00  0.9888528E-01  0.4721481E-01
#k= 2 12.00 days  l1,l2,linf=  0.2320193E+01  0.1051828E+01  0.4538982E+00
#k= 3 12.00 days  l1,l2,linf=  0.4281969E+01  0.1469139E+01  0.1089794E+01
#k= 4 12.00 days  l1,l2,linf=  0.1001085E+01  0.9986742E+00  0.9999014E+00
#
#
###############################################################
# Leapfrog default benchmark (used to check nothing is broken)
###############################################################
# NE=11, dt=360, nu=0  limiter=0  filter_freq=0 NV=4
# set smooth = 0.05
# set integration = explicit
#
# Error for cosine bell & other solid body rotation problems:
#k= 1 12.00 days  l1,l2,linf=  0.3080480E+00  0.9562703E-01  0.4883010E-01
#k= 2 12.00 days  l1,l2,linf=  0.2319669E+01  0.1049853E+01  0.4539350E+00
#k= 3 12.00 days  l1,l2,linf=  0.3653317E+01  0.1378439E+01  0.1020411E+01
#k= 4 12.00 days  l1,l2,linf=  0.1000506E+01  0.9986733E+00  0.9998071E+00
#
# 
set wdir = ~/scratch1/swtc1
set src = ~/codes/homme/build.Linux
set input = ~/codes/homme/test/sw_conservative
set NCPU = 16

if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
if ( ${?SLURM_NNODES} ) then
   # SLURM_NNODES  = number of nodes
   # hard to tell how many cores per nodes
   # set NCPU to zero, and mpirun will use the max allowed
   set NCPU = 0
endif
echo NCPU = $NCPU


set test_case = swtc1
set params = $input/Params4.inc      # 4 levels
set SUBCASE = 
set ndays = 12

diff  $params  $src/../Params.inc
if ($status != 0) then
   echo "replacing Params.inc"
   cp $params $src/../Params.inc
endif



cd $src
rm -f sweqx
make -j2 sweqx
mkdir $wdir
cd $wdir
mkdir movies


set NE=11
set tstep = 360
set nu = 0
set filter_freq = 0

set integration = explicit
#set integration = runge_kutta

if ( $integration == 'explicit') then
   # LEAPFROG SETTINGS - to mimic dynamics used in 3D
   set smooth = 0.05 ; set LFTfreq = 0   # pure LF  used in HOMME 3D w/o subcycling
   #set smooth = 0 ; set LFTfreq = 4       # RK2-4  used in HOMME 3D subcyclign
   set rk_stage = 0
   set limiter = 0  # leapfrog code cannot use limiters
endif

if ( $integration == 'runge_kutta') then
   # RKSSP settings - to mimic tracer advection used in 3D
   set smooth = 0
   set LFTfreq = 0
   set rk_stage = 3
   set limiter = 0
endif




set name = ${test_case}-NE${NE}-t${tstep}-limiter$limiter

set sfreq = 24
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


set nu_s = $nu
sed s/^ne.\*/"ne = $NE"/  $input/swtc1.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
sed s/limiter_option.\*/"limiter_option = $limiter"/  |\
sed s/smooth.\*/"smooth = $smooth"/  |\
sed s/test_case.\*/"test_case = \'$test_case\'   $SUBCASE"/  |\
sed s/explicit/$integration/  |\
sed s/rk_stage_user.\*/"rk_stage_user = $rk_stage"/  |\
sed s/LFTfreq.\*/"LFTfreq = $LFTfreq"/  |\
sed s/ndays.\*/"ndays = $ndays"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu_s"/  |\
sed s/filter_freq.\*/"filter_freq = $filter_freq"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> swtc1.nl

date
mpirun -np $NCPU $src/sweqx < swtc1.nl | tee  sweq.out
date

