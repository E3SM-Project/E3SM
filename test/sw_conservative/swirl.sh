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
set OUTUNITS = 1  
set OUTFREQ =  0

set test_case = swirl
set ndays = 5
set SUBCASE = "sub_case = 4"

set build = 0
set make = 0
if ( $build == 1 ) then
   cd $src
   ./configure --enable-blas --enable-lapack --with-netcdf=$NETCDF_PATH \
    --with-pnetcdf=$PNETCDF_PATH NP=4 PLEV=4   --enable-energy-diagnostics
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
set limiter = 8


# convergence study 0,8,81

#
# limiter_option=0
# limiter_option=4     original 0 limiter
# limiter_option=8     prevent undershoot and overshoot
# limiter_option=81    prevent undershoot
# limiter_option=84    minimum 0/0.1 hard coded
#

set NE = 160

if ( $NE == 10 ) then
   set nu = 256e14  # dx^3 scaled from NE=40 4e14 ref value
   set tstep = 360  # default
   #set tstep = 180  # default
endif

if ( $NE == 20) then
   #set nu = 6.6e14    # ref value, used for dx^4 scaling
   set nu = 32e14  # dx^3 scaled from NE=40 4e14 ref value    
   set tstep = 180  # default

   # slotted cylinder, lim8  day 5 min/max   0.999999999700422E-01   0.996964667304424E+00 
# 1400: ok
# 1425: ok
endif


if ( $NE == 27) then
   # NE=26 errors got as low as 0.0344.  
   # for N&L equiv. resolution, had to go to NE=27 for < 0.033
   # lim81, tstep=120, nu=6e13  l2error = 0.031   tstep=60  l2error=.031
   # lim81, tstep=120, nu=9e13  l2error = 0.0309
   # lim81, tstep=120, nu=1e14  l2error = 0.0310
   # lim81, tstep=120, nu=2e14  l2error = 0.0326
   # lim81, tstep=120, nu=3e14  l2error = 0.0356
   set nu = 2e14    # scaled from NE=20 6.6e14
   set tstep = 120  #default
endif

if ( $NE == 30 ) then
  # Test case describes CFL = dt Umax / dx    
  #   dt(days), Umax=3.26 radians/day, dx(radians)   CFL=0.26
  # OR:  Rearth=6.376e6  
  #   dt(seconds) * 3.26*Rearth/(24*3600) / ( 2*pi*Rearth * deg/360)
  #   dt(seconds) 240.6(m/s) / 111,282m :: CFL=0.26
  #
  # CAM dynamics:  tstep=90  based on 340m/s gravity wave   CFL=0.27 (leapfrog)
  # swirl test case winds:  240m/s.  RKSSP 3 stage has sqrt(2.0) larger CFL: dt=180
  #  
  set tstep = 120   # CFL=0.26
  set nu = 1.3e14    # scaled from NE=20 6.6e14

#  set tstep = 133.3333   # CFL = .288   SSP at day 5. 
#  set tstep = 144   # CFL=0.311   going unstable (1e2) at day 5
#  set tstep = 150   # CFL=0.324   NaN  (at day 5)

# new scheme
#   set tstep=150    #       l2error=.046
#   set tstep=240    #       l2error=.060
#   set tstep=300    #       l2error=.070
#   set tstep=360    #       l2error=.060   NAN
#   set nu = 0
endif


if ( $NE == 40) then
   # lim81, tstep=90, nu=1.25e13:   l2error=.0137        
   # lim81, tstep=90, nu=4e13:      l2error=.0114
   # lim81, tstep=90, nu=8e13:      l2error=.0108    tstep=45:  0.1089802E-01
   # lim81, tstep=90, nu=1e14:      l2error=.0111
   # lim81, tstep=90, nu=5e14:      l2error=.0283
   #set nu = 4e13     # scaled from NE=20 6.6e14
   set nu = 4e14  # dx^3 scaled from NE=40 4e14 ref value
   set tstep = 90  # default 
   #set tstep = 45
endif
if ( $NE == 80) then
  # grid spacing:  41,730m   Umax=240.6m/s
   #set tstep = 45  # unstable        CFL=.26
   #set tstep = 40  # errors are bad  CFL=.23
   set tstep = 30   # new default    CFL=.17
   #set tstep = 15   
   #set tstep = 3
   #set nu = 2.5e12    # scaled from NE=20 6.6e14
   set nu = 5e13  # dx^3 scaled from NE=40 4e14 ref value
endif
if ( $NE == 160) then
   #set tstep = 15
   set tstep = 6
   #set tstep = 2
   #set nu = 1.6e11  # scaled from NE=20 6.6e14
   set nu = 6.25e12  # dx^3 scaled from NE=40 4e14 ref value
endif




set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}-nu$nu-limiter$limiter


set sfreq = 12
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




