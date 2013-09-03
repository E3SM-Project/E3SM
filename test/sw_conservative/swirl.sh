#!/bin/tcsh -f
#
# 50 nodes, 1h, NE=80  dt=3
# 50 nodes, 2h, NE=160 dt=5  
# 100 nodes, 4h, NE=160 dt=1  
#
#SBATCH --job-name swirl
#SBATCH -N 40
#SBATCH --account=FY139209
#SBATCH --time=1:00:00
#


set wdir = ~/scratch1/swirl
set src = ~/codes/homme/build/sweqx
set input = ~/codes/homme/test/sw_conservative
set NCPU = 1

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


set build = 0
set make = 0
set NP = 7
set limiter = 0


set nu = 0
if ( $NP == 4 ) then
# convergence study 0,8,81
   # NE=26 errors got as low as 0.0344.  
   # for N&L equiv. resolution, had to go to NE=27 for < 0.033
   # lim81, tstep=120, nu=6e13  l2error = 0.031   tstep=60  l2error=.031
   # lim81, tstep=120, nu=9e13  l2error = 0.0309
   # lim81, tstep=120, nu=1e14  l2error = 0.0310
   # lim81, tstep=120, nu=2e14  l2error = 0.0326
   # lim81, tstep=120, nu=3e14  l2error = 0.0356
   #set nu = 2e14    # scaled from NE=20 6.6e14
   #set tstep = 120  #default

   if ( $limiter == 0 ) then
     set NE = 10 ; set tstep = 360   ; set nu=1.1e16
     #set NE = 20 ; set tstep =  90   ; set nu=6.6e14
     #set NE = 40 ; set tstep =  22.5 ; set nu=4.1e13
     #set NE = 80 ; set tstep =  6    ; set nu=2.6e12
   endif
   if ( $limiter == 84 ) then
     set NE = 10 ; set tstep = 360   ; set nu=1.1e16
     #set NE = 20 ; set tstep =  90   ; set nu=6.6e14
     #set NE = 40 ; set tstep =  22.5 ; set nu=4.1e13
     #set NE = 80 ; set tstep =  6    ; set nu=2.6e12
   endif
   if ( $limiter == 8 ) then
     set NE = 10 ; set tstep = 360   ; set nu=1.1e16
     #set NE = 20 ; set tstep =  120   ; set nu=6.6e14
     #set NE = 40 ; set tstep =  45 ; set nu=4.1e13
     #set NE = 80 ; set tstep =  16    ; set nu=2.6e12
   endif
endif

if ( $NP == 7 ) then
  if ( $limiter == 0) then
    set NE = 5  ; set tstep = 270 ; set nu=1.30e16
    #set NE = 5  ; set tstep = 150 ; set nu=1.30e16
    #set NE = 10 ; set tstep = 75 ; set nu=1.00e14
    #set NE = 20 ; set tstep =  10 ; set nu=7.80e11   #l2: cos .57e-2  gauss: .6e-4
    #set NE = 20 ; set tstep =  5 ; set nu=7.80e11     #l2: cos .57e-2  gauss:  .58e-4
    #set NE = 40 ; set tstep =  1 ; set nu=6.10e09
  endif
  if ( $limiter == 8) then
    #set NE = 5  ; set tstep = 270 ; set nu=1.30e16
    #set NE = 10 ; set tstep = 100 ; set nu=1.00e14
    #set NE = 20 ; set tstep =  30 ; set nu=7.80e11
    set NE = 40 ; set tstep =  10 ; set nu=6.10e09
  endif
  if ( $limiter == 84 ) then
    set NE = 5  ; set tstep = 270 ; set nu=1.30e16
    #set NE = 10 ; set tstep = 67.5 ; set nu=1.00e14
    #set NE = 20 ; set tstep =  15 ; set nu=7.80e11
    #set NE = 40 ; set tstep =  4 ; set nu=6.10e09
  endif
endif



if ( $build == 1 ) then
   cd $src
   ./configure   NP=$NP PLEV=4  \
    --with-netcdf=$NETCDF_PATH --with-pnetcdf=$PNETCDF_PATH
   make depends
   make -j4 clean
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



# output units: 0,1,2 = timesteps, days, hours
set OUTUNITS = 2  
set OUTFREQ =  30  # output 0,1.25,2.5,3.75,5.0 days

set test_case = swirl
set ndays = 5
set SUBCASE = "sub_case = 4"


# RKSSP settings - to mimic tracer advection used in 3D
set smooth = 0
set integration = runge_kutta
set LFTfreq = 0
set rk_stage = 3
set filter_freq = 0
set sfreq = 12
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


set name = ne${NE}np${NP}-t${tstep}-nu$nu-lim$limiter
set wdir = $wdir/$name
mkdir $wdir
cd $wdir
mkdir movies




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
sed s/output_frequency.\*/"output_frequency = $OUTFREQ"/  |\
sed s/output_timeunits.\*/"output_timeunits = $OUTUNITS  interp_type=1"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> input.nl

date
mpirun -np $NCPU $src/sweqx.swirl < input.nl | tee  sweq.out
date
# if timing was turned on, get sweq cost:
grep sweq HommeSWTime | sort | tail -1

#  ncl $input/geo.ncl




