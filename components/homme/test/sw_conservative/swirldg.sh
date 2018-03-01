#!/bin/tcsh -f
#
# 50 nodes, 1h, NE=80  dt=3
# 50 nodes, 2h, NE=160 dt=5  
# 100 nodes, 4h, NE=160 dt=1  
#
#SBATCH --job-name swirl
#SBATCH -N 60
#SBATCH --account=FY127788
#SBATCH --time=2:00:00
#XXXTCH --depend=afterany:jobid
#
#XPBS -l nodes=200:ppn=2
#XPBS -l walltime=8:00:00
#XPBS -N swirl
#XPBS -j oe
#XPBS -A FY081407
#XPBX -a 1758
#XXX -W depend=afterany:jobid


set wdir = ~/scratch1/swirl
set src = ~/codes/homme/build/sweqx
set input = ~/codes/homme/test/sw_conservative
set NCPU = 1
set NP = 4

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

# output units: 0,1,2 = timesteps, days, hours
set OUTUNITS = 2  
set OUTFREQ =  4

set test_case = swirl
set ndays = 5
set SUBCASE = "sub_case = 4"

set build = 0
set make = 1
if ( $build == 1 ) then
   cd $src
   ./configure --enable-blas --enable-lapack --with-netcdf=$NETCDF_PATH \
    --with-pnetcdf=$PNETCDF_PATH NP=$NP PLEV=4   --enable-energy-diagnostics
   make depends
   make clean
   make -j4 sweqx
   mv sweqx sweqx.swirl
   exit
endif
if ( $make == 1 ) then
   cd src
   make -j4 sweqx
   if ($status) exit
   mv sweqx sweqx.swirl
endif


mkdir $wdir
cd $wdir
mkdir movies

set nu = 0


# LEAPFROG SETTINGS - to mimic dynamics used in 3D
#set smooth = 0.05 ; set LFTfreq = 0   # pure LF  used in HOMME 3D w/o subcycling
#set smooth = 0 ; set LFTfreq = 4       # RK2-4  used in HOMME 3D subcyclign
#set rk_stage = 0
#set integration = explicit
#set limiter = 0


# RKSSP settings - to mimic tracer advection used in 3D
set smooth = 0
set integration = runge_kutta
set LFTfreq = 0
set rk_stage = 3
set limiter = 4


# convergence study 0,8,81

#
# limiter_option=0
# limiter_option=4     original 0 limiter
# limiter_option=8     prevent undershoot and overshoot
# limiter_option=81    prevent undershoot
# limiter_option=84    minimum 0/0.1 hard coded
#

set NE = 30
set tstep=120


if ( $NE == 10 ) then
   set tstep = 360  # default
   #set tstep = 180  # default
endif

if ( $NE == 20) then
   set tstep = 180  # default
endif


if ( $NE == 27) then
   # NE=26 errors got as low as 0.0344.  
   # for N&L equiv. resolution, had to go to NE=27 for < 0.033
   # lim81, tstep=120, nu=6e13  l2error = 0.031   tstep=60  l2error=.031
   # lim81, tstep=120, nu=9e13  l2error = 0.0309
   # lim81, tstep=120, nu=1e14  l2error = 0.0310
   # lim81, tstep=120, nu=2e14  l2error = 0.0326
   # lim81, tstep=120, nu=3e14  l2error = 0.0356
   set tstep = 120  #default
endif

if ( $NE == 30 ) then
  set tstep = 120   # CFL=0.26
endif


if ( $NE == 40) then
   set tstep = 90  # default 
endif
if ( $NE == 80) then
   set tstep = 30   # new default    CFL=.17
endif
if ( $NE == 160) then
   set tstep = 6
endif




set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}-limiter$limiter


set sfreq = 4
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


set nu_s = $nu
sed s/^ne.\*/"ne = $NE  npdg=4"/  $input/swtc1.nl |\
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
sed s/output_frequency.\*/"output_frequency = $OUTFREQ"/  |\
sed s/output_timeunits.\*/"output_timeunits = $OUTUNITS"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> input.nl

date
mpirun -np $NCPU $src/sweqx.swirl < input.nl | tee  sweq.out
date
# if timing was turned on, get sweq cost:
grep sweq HommeSWTime | sort | tail -1

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out
mv -f swirl.l1.errors $name.l1
mv -f swirl.l2.errors $name.l2
mv -f swirl.linf.errors $name.linf
mv -f movies/swirl1.nc movies/$name.nc
#  ncl $input/geo.ncl




