#!/bin/tcsh -f
#
# 50 nodes, 1h, NE=80  dt=3
# 50 nodes, 2h, NE=160 dt=5  
# 100 nodes, 4h, NE=160 dt=1  
#
#SBATCH --job-name swirl
#SBATCH -N 4
#SBATCH --account=FY150001
#SBATCH -p ec
#SBATCH --time=0:10:00
#

set wdir = ~/scratch1/sweqx
set HOMME = ~/codes/homme
#set MACH = $HOMME/cmake/machineFiles/rhel5.cmake
set MACH = $HOMME/cmake/machineFiles/redsky.cmake
set input = $HOMME/test/sw_conservative

set builddir = $wdir/bld
set rundir = $wdir/swtc2ref
mkdir -p $rundir
mkdir -p $wdir/bld
cd $builddir

set NP = 4
set limiter = 8


set NCPU = 4
if ( ${?SLURM_NNODES} ) then
   # set NCPU to zero, and mpirun will use the max allowed
   set NCPU = 0
endif
echo NCPU = $NCPU


#configure the model
set conf = 0
set make = 1
#cmake:
cd $builddir
if ( $conf == 1 ) then
   rm -rf CMakeFiles CMakeCache.txt
   cmake -C $MACH -DSWEQX_PLEV=6  -DSWEQX_NP=$NP $HOMME
   make -j4 clean
   make -j4 sweqx
   exit
endif
if ( $make == 1 ) then
   make -j4 sweqx
    if ($status) exit
endif
set exe = $builddir/src/sweqx/sweqx



cd $rundir
mkdir movies


# NOTE:
# shallow water test was converted to use new/physical units 7/2015, and the
# case now runs for 12 days (instead of 5 days).  
# 
# settings below given in 'OLD' units need to be converted.
#   tstep should be INCREASED by 12/5
#   viscosity should be DECREASED by 12/5
# 


set nu = 0
if ( $NP == 4 ) then
# convergence study 0,8,81
   # NE=26 errors got as low as 0.0344.  
   # for N&L equiv. resolution, had to go to NE=27 for < 0.033
   # OLD UNITS:
   # lim81, tstep=120, nu=6e13  l2error = 0.031   tstep=60  l2error=.031
   # lim81, tstep=120, nu=9e13  l2error = 0.0309
   # lim81, tstep=120, nu=1e14  l2error = 0.0310
   # lim81, tstep=120, nu=2e14  l2error = 0.0326
   # lim81, tstep=120, nu=3e14  l2error = 0.0356

   # settings below have been converted to new/physical units:
   if ( $limiter == 0 ) then
     # OLD units: NEEDS UPDATE: 
     #set NE = 10 ; set tstep = 360   ; set nu=1.1e16
     #set NE = 20 ; set tstep =  90   ; set nu=6.6e14
     #set NE = 40 ; set tstep =  22.5 ; set nu=4.1e13
     #set NE = 80 ; set tstep =  6    ; set nu=2.6e12
   endif
   if ( $limiter == 84 ) then
     # OLD units: NEEDS UPDATE: 
     #set NE = 10 ; set tstep = 360   ; set nu=1.1e16
     #set NE = 20 ; set tstep =  90   ; set nu=6.6e14
     #set NE = 40 ; set tstep =  22.5 ; set nu=4.1e13
     #set NE = 80 ; set tstep =  6    ; set nu=2.6e12
   endif
   if ( $limiter == 8 ) then
     # OLD UNITS, DONT USE
     #set NE = 10 ; set tstep = 360   ; set nu=1.1e16   #l2erros:  .349  .183 .437
     #set NE = 20 ; set tstep =  120   ; set nu=6.6e14  #l2errors: .121  .0461  .331
     #set NE = 40 ; set tstep =  45 ; set nu=4.1e13
     #set NE = 80 ; set tstep =  16    ; set nu=2.6e12

     # NEW units, scaled to match above:
     #set NE = 10 ; set tstep = 720   ; set nu=4.6e15   #l2errors: .347 .182  .437
     set NE = 20 ; set tstep = 240   ; set nu=2.8e14    #l2errors: .120 .0458  .331 
     #set NE = 40 ; set tstep = 100  ; set nu=1.7e13  
     #set NE = 80 ; set tstep = 36   ; set nu=1.1e12  
   endif
endif

if ( $NP == 7 ) then
  if ( $limiter == 0) then
    # OLD units: NEEDS UPDATE: 
    #set NE = 5  ; set tstep = 270 ; set nu=1.30e16
    #set NE = 5  ; set tstep = 150 ; set nu=1.30e16
    #set NE = 10 ; set tstep = 75 ; set nu=1.00e14
    #set NE = 20 ; set tstep =  10 ; set nu=7.80e11   #l2: cos .57e-2  gauss: .6e-4
    #set NE = 20 ; set tstep =  5 ; set nu=7.80e11     #l2: cos .57e-2  gauss:  .58e-4
    #set NE = 40 ; set tstep =  1 ; set nu=6.10e09
  endif
  if ( $limiter == 8) then
    # OLD units: NEEDS UPDATE: 
    #set NE = 5  ; set tstep = 270 ; set nu=1.30e16
    #set NE = 10 ; set tstep = 100 ; set nu=1.00e14
    #set NE = 20 ; set tstep =  30 ; set nu=7.80e11
    #set NE = 40 ; set tstep =  10 ; set nu=6.10e09
  endif
  if ( $limiter == 84 ) then
    # OLD units: NEEDS UPDATE: 
    #set NE = 5  ; set tstep = 270 ; set nu=1.30e16
    #set NE = 10 ; set tstep = 67.5 ; set nu=1.00e14
    #set NE = 20 ; set tstep =  15 ; set nu=7.80e11
    #set NE = 40 ; set tstep =  4 ; set nu=6.10e09
  endif
endif



# output units: 0,1,2 = timesteps, days, hours
set OUTUNITS = 2  
set OUTFREQ =  30  # output 0,1.25,2.5,3.75,5.0 days

set test_case = swirl
set ndays = 12
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
mpirun -np $NCPU $exe < input.nl | tee  sweq.out
date
# if timing was turned on, get sweq cost:
grep sweq HommeSWTime | sort | tail -1

#  ncl $input/geo.ncl




