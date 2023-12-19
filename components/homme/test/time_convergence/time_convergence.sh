#!/bin/tcsh 
#
#   Jobscript for launching dcmip2012 test3-1 on SNL clusters
#SBATCH -J d31-theta          # job name              
#SBATCH -N 1                  # total number of mpi tasks requested 
#SBATCH -p ec                 # queue (partition) -- normal, development, etc. 
#SBATCH -t 00:10:00           # run time (hh:mm:ss)   
#SBATCH --account=FY210162                                           

set refsol_namelist  = namelist-refsol.nl
set testsol_namelist = namelist-testsol.nl

echo message | ts

set OMP_NUM_THREADS = 1
set NCPU = 1

if ( ${?PBS_ENVIRONMENT} ) then   # anvil      
  set NCPU = $PBS_NNODES
  if ( $PBS_ENVIRONMENT == PBS_BATCH ) cd $PBS_O_WORKDIR
endif
if ( ${?SLURM_NNODES} ) then   # skybridge                
    set NCPU = $SLURM_NNODES
    @ NCPU *= 16
    @ NCPU /= $OMP_NUM_THREADS
endif

#  set paths to source code, build directory and run directory #                                   
set blddir  =  ~/e3sm_work/test_convergence # directory in which you will build homme
set testdir = $PWD # input files for test case, I assume you start in the test_convergence directory in the tests                    
set HOMME   = $PWD/../..  # HOMME base directory is two levels up from the time_convergence testdir

set MACH = $HOMME/cmake/machineFiles/sandia-srn-sems.cmake            
#set MACH = $HOMME/cmake/machineFiles/anvil.cmake
#set MACH = $HOMME/cmake/machineFiles/titan.cmake                                                   
#set MACH = $HOMME/cmake/machineFiles/titan-cuda.cmake    # CMAKE file for Titan GPU support                
#set MACH = $HOMME/cmake/machineFiles/darwin.cmake      

# vertical resolution is 20 in dcmip3.1
set nlev  = 20       # vertical resolution                                                                           
set qsize = 0     # number of passive tracers     

cd $blddir
echo $PWD
set exe = $blddir
set build = 1  # set to 1 to force build                                                                     
# rm $bld/CMakeCache.txt    # remove this file to force re-configure                                           
if (! -f CMakeCache.txt) then
   rm -rf CMakeFiles CMakeCache.txt
   echo "running CMAKE to configure the model"
   cmake -C $MACH \
   -DQSIZE_D=$qsize \
   -DPREQX_PLEV=$nlev \
   -DPREQX_NP=4  \
   -DPREQX_USE_PIO=TRUE                          \
   -DPREQX_USE_ENERGY=TRUE \
   -DBUILD_HOMME_PREQX_KOKKOS="OFF" \
   -DBUILD_HOMME_PREQX_ACC="OFF" \
   -DHOMME_ENABLE_COMPOSE="OFF" \
   -DPREQX_USE_ENERGY=TRUE \
    $HOMME
#   make -j4 clean
endif
if ( ! -f $exe) set build = 1   # no exe, force build                           
if ( $build == 1 ) then
   make -j4 $blddir
   if ($status) exit
endif

make -j32 -f Makefile theta-l-nlev20-native

cd dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta-l

cp $testdir/$refsol_namelist .
cp $testdir/$testsol_namelist .

set EXEC = $blddir/test_execs/theta-l-nlev20-native/theta-l-nlev20-native

#############################################################################
# Construct reference solution
############################################################################# 
mpirun -np $NCPU $EXEC < $testdir/$refsol_namelist


######################################
#mpirun -np $NCPU $EXEC < input-testsol.nl



echo message | ts
