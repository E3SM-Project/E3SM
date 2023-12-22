#!/bin/tcsh                                                                                                                                                        
#                                                                                                                                                                  
#   Jobscript for launching dcmip2012 test3-1 on SNL clusters                                                                                                      
#SBATCH -J d31-theta          # job name                                                                                                                           
#SBATCH -N 2                  # total number of mpi tasks requested                                                                                               
#SBATCH -p ec                 # queue (partition) -- normal, development, etc.                                                                                     
#SBATCH -t 01:00:00           # run time (hh:mm:ss)                                                                                                                
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
set blddir    = '~/e3sm_work/test_convergence' # directory in which you will build homme
set testdir   = $PWD # input files for test case, I assume you start in the test_convergence directory in the tests                    
set HOMME     = $PWD/../..  # HOMME base directory is two levels up from the time_convergence testdir
set outputdir = '/nscratch/asteyer/test_convergence/' # this is where the output will go 


set MACH = $HOMME/cmake/machineFiles/sandia-srn-sems.cmake            
#set MACH = $HOMME/cmake/machineFiles/anvil.cmake
#set MACH = $HOMME/cmake/machineFiles/titan.cmake                                                   
#set MACH = $HOMME/cmake/machineFiles/titan-cuda.cmake    # CMAKE file for Titan GPU support                
#set MACH = $HOMME/cmake/machineFiles/darwin.cmake      

# vertical resolution is 20 in dcmip3.1
set nlev  = 20       # vertical resolution                                                                           
set qsize = 0     # number of passive tracers     

cd $blddir
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
   -DBUILD_HOMME_SWDGX=FALSE                     \
   -DBUILD_HOMME_SWEQX=FALSE                     \
   -DBUILD_HOMME_PRIMDGX=FALSE                   \
   -DPREQX_USE_PIO=TRUE                          \
   -DPREQX_USE_ENERGY=TRUE \
   -DBUILD_HOMME_PREQX_KOKKOS="OFF" \
   -DBUILD_HOMME_PREQX_ACC="OFF" \
   -DHOMME_ENABLE_COMPOSE="OFF" \
    $HOMME
#   make -j4 clean
endif
if ( ! -f $exe) set build = 1   # no exe, force build                           
if ($build == 1 ) then
   make -j4 $blddir
   if ($status) exit
endif

make -j32 -f Makefile theta-l-nlev20-native

cd $blddir/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta-l
make -j32 install

cp $testdir/$refsol_namelist .
cp $testdir/$testsol_namelist .

set EXEC = $blddir/test_execs/theta-l-nlev20-native/theta-l-nlev20-native

#############################################################################                                                                                      
#############################################################################                                                                                      
#############################################################################
# Construct nonhydrostatic reference solution
############################################################################# 
#############################################################################                                                                                      
#############################################################################                                                                                      
set stroutput = "output_dir"
cp $testdir/$refsol_namelist refnl.nl
set linenum = `grep -n $stroutput refnl.nl | cut -d: -f1`
set new = '  output_dir ="'$outputdir'"'
sed -i "${linenum}s,.*,$new," refnl.nl
echo 'Computing nonhydrostatic reference solution with dt = .005 ...'
mpirun -np $NCPU $EXEC < refnl.nl 
mv $outputdir/dcmip2012_test31.nc $outputdir/refnhsol.nc
echo 'Done computing nonhydrostatic reference solution'

set stroutput = "output_dir"
cp $testdir/$refsol_namelist refnl.nl
set linenum = `grep -n $stroutput refnl.nl | cut -d: -f1`
set new = '  output_dir ="'$outputdir'"'
sed -i "${linenum}s,.*,$new," refnl.nl


set stroutput = "output_dir"
echo $stroutput
cp $testdir/$testsol_namelist testnl.nl
set linenum = `grep -n $stroutput testnl.nl | cut -d: -f1`
set new = '  output_dir ="'$outputdir'"'                          
sed -i "${linenum}s,.*,$new," testnl.nl  

#############################################################################                                                                                      
#############################################################################                                                                                      
##### tstep_type = 9 NH RUNS
#############################################################################                                                                                      
#############################################################################    
###########################################################                                                                
##### tstep_type = 9, dt = 3.6, nonhydrostatic        #####                                                                                     
###########################################################    
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=8"
set new3 = "  output_frequency=8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=3.6"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl                                                                                                
echo 'Computing tstep_type = 9, dt = 3.6...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9_dt3p6.nc
echo 'Done computing tstep_type = 9, dt = 3.6...'

###########################################################                                                                                                        
##### tstep_type = 9, dt = 1.8, nonhydrostatic        #####                                                                                                        
###########################################################                                                                                                  
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=16"
set new3 = "  output_frequency=16"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=1.8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 9, dt = 1.8...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9_dt1p8.nc
echo 'Done computing tstep_type = 9, dt = 1.8...'

###########################################################                                                                                                        
##### tstep_type = 9, dt = 0.9, nonhydrostatic        #####                                                                                                        
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=32"
set new3 = "  output_frequency=32"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=0.9"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 9, dt = 0.9...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9_dt0p9.nc
echo 'Done computing tstep_type = 9, dt = 0.9...'

###########################################################                                                                                                        
##### tstep_type = 9, dt = 0.45, nonhydrostatic       #####                                                                                                       
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=64"
set new3 = "  output_frequency=64"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=0.45"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 9, dt = 0.45...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9_dt0p45.nc
echo 'Done computing tstep_type = 9, dt = 0.45...'

#############################################################################                                                                                      
#############################################################################                                                                                      
##### tstep_type = 10 NH RUNS                                                                                                                    
#############################################################################                                                                                      
#############################################################################    
###########################################################                                                                
##### tstep_type = 10, dt = 3.6, nonhydrostatic       #####                                                                                     
###########################################################    
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=8"
set new3 = "  output_frequency=8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=3.6"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl                                                                                                
echo 'Computing tstep_type = 10, dt = 3.6...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10_dt3p6.nc
echo 'Done computing tstep_type = 10, dt = 3.6...'

###########################################################                                                                                                        
##### tstep_type = 10, dt = 1.8, nonhydrostatic       #####                                                                                                        
###########################################################                                                                                                  
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=16"
set new3 = "  output_frequency=16"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=1.8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 10, dt = 1.8...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10_dt1p8.nc
echo 'Done computing tstep_type = 10, dt = 1.8...'

###########################################################                                                                                                        
##### tstep_type = 10, dt = 0.9, nonhydrostatic       #####                                                                                                        
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=32"
set new3 = "  output_frequency=32"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=0.9"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 10, dt = 0.9...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10_dt0p9.nc
echo 'Done computing tstep_type = 10, dt = 0.9...'

###########################################################                                                                                                        
##### tstep_type = 10, dt = 0.45, nonhydrostatic      #####                                                                                                       
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=64"
set new3 = "  output_frequency=64"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=0.45"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 10, dt = 0.45...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10_dt0p45.nc
echo 'Done computing tstep_type = 10, dt = 0.45...'

#############################################################################                                                                                      
#############################################################################                                                                                      
##### tstep_type = 5 NH RUNS                                                                                                                                 
#############################################################################                                                                                      
############################################################################# 
###########################################################                                                                
##### tstep_type = 5, dt = 3.6, nonhydrostatic        #####                                                                                     
###########################################################    
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=8"
set new3 = "  output_frequency=8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=3.6"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl                                                                                                
echo 'Computing tstep_type = 5, dt = 3.6...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5_dt3p6.nc
echo 'Done computing tstep_type = 5, dt = 3.6...'

###########################################################                                                                                                        
##### tstep_type = 5, dt = 1.8, nonhydrostatic        #####                                                                                                        
###########################################################                                                                                                  
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=16"
set new3 = "  output_frequency=16"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=1.8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 5, dt = 1.8...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5_dt1p8.nc
echo 'Done computing tstep_type = 5, dt = 1.8...'

###########################################################                                                                                                        
##### tstep_type = 5, dt = 0.9, nonhydrostatic        #####                                                                                                        
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=32"
set new3 = "  output_frequency=32"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=0.9"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 5, dt = 0.9...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5_dt0p9.nc
echo 'Done computing tstep_type = 5, dt = 0.9...'

###########################################################                                                                                                        
##### tstep_type = 5, dt = 0.45, nonhydrostatic       #####                                                                                                       
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=64"
set new3 = "  output_frequency=64"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=0.45"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 5, dt = 0.45...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5_dt0p45.nc
echo 'Done computing tstep_type = 5, dt = 0.45...'


#############################################################################                                                                                      
#############################################################################                                                                                      
#############################################################################                                                                                      
# Construct hydrostatic reference solution                                                                                                                      
#############################################################################                                                                                      
#############################################################################                                                                                      
#############################################################################     
set strh = "hydrostatic_mode"
set linenum = `grep -n $strh refnl.nl | cut -d: -f1`
set new = '  theta_hydrostatic_mode =.true.'
sed -i "${linenum}s,.*,$new," refnl.nl
echo 'Computing hydrostatic reference solution with dt = .005 ...'
mpirun -np $NCPU $EXEC < refnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/refhsol.nc
echo 'Done computing hydrostatic reference solution'


set strh = "hydrostatic_mode"
set linenum = `grep -n $strh testnl.nl | cut -d: -f1`
set new = '  theta_hydrostatic_mode =.true.'
sed -i "${linenum}s,.*,$new," testnl.nl

###########################################################                                                                
##### tstep_type = 9, dt = 3.6, hydrostatic           #####                                                                                     
###########################################################    
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=8"
set new3 = "  output_frequency=8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=3.6"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl                                                                                                
echo 'Computing tstep_type = 9, dt = 3.6...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9h_dt3p6.nc
echo 'Done computing tstep_type = 9, dt = 3.6...'

###########################################################                                                                                                        
##### tstep_type = 9, dt = 1.8, hydrostatic           #####                                                                                                        
###########################################################                                                                                                  
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=16"
set new3 = "  output_frequency=16"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=1.8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 9, dt = 1.8...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9h_dt1p8.nc
echo 'Done computing tstep_type = 9, dt = 1.8...'

###########################################################                                                                                                        
##### tstep_type = 9, dt = 0.9, hydrostatic           #####                                                                                                        
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=32"
set new3 = "  output_frequency=32"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=0.9"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 9, dt = 0.9...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9h_dt0p9.nc
echo 'Done computing tstep_type = 9, dt = 0.9...'

###########################################################                                                                                                        
##### tstep_type = 9, dt = 0.45, hydrostatic          #####                                                                                                       
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=9"
set new2 = "  nmax=64"
set new3 = "  output_frequency=64"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=9"
set new2 = "  tstep=0.45"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 9, dt = 0.45...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep9h_dt0p45.nc
echo 'Done computing tstep_type = 9, dt = 0.45...'

#############################################################################                                                                                      
#############################################################################                                                                                      
##### tstep_type = 10 NH RUNS                                                                                                                    
#############################################################################                                                                                      
#############################################################################    
###########################################################                                                                
##### tstep_type = 10, dt = 3.6, hydrostatic          #####                                                                                     
###########################################################    
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=8"
set new3 = "  output_frequency=8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=3.6"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl                                                                                                
echo 'Computing tstep_type = 10, dt = 3.6...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10h_dt3p6.nc
echo 'Done computing tstep_type = 10, dt = 3.6...'

###########################################################                                                                                                        
##### tstep_type = 10, dt = 1.8, hydrostatic          #####                                                                                                        
###########################################################                                                                                                  
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=16"
set new3 = "  output_frequency=16"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=1.8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 10, dt = 1.8...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10h_dt1p8.nc
echo 'Done computing tstep_type = 10, dt = 1.8...'

###########################################################                                                                                                        
##### tstep_type = 10, dt = 0.9, hydrostatic          #####                                                                                                        
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=32"
set new3 = "  output_frequency=32"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=0.9"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 10, dt = 0.9...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10h_dt0p9.nc
echo 'Done computing tstep_type = 10, dt = 0.9...'

###########################################################                                                                                                        
##### tstep_type = 10, dt = 0.45, hydrostatic         #####                                                                                                       
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=10"
set new2 = "  nmax=64"
set new3 = "  output_frequency=64"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=10"
set new2 = "  tstep=0.45"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 10, dt = 0.45...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep10h_dt0p45.nc
echo 'Done computing tstep_type = 10, dt = 0.45...'


#############################################################################                                                                                      
#############################################################################                                                                                      
##### tstep_type = 5 NH RUNS                                                                                                                                 
#############################################################################                                                                                      
############################################################################# 
###########################################################                                                                
##### tstep_type = 5, dt = 3.6, hydrostatic           #####                                                                                     
###########################################################    
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=8"
set new3 = "  output_frequency=8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=3.6"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl                                                                                                
echo 'Computing tstep_type = 5, dt = 3.6...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5h_dt3p6.nc
echo 'Done computing tstep_type = 5, dt = 3.6...'

###########################################################                                                                                                        
##### tstep_type = 5, dt = 1.8, hydrostatic           #####                                                                                                        
###########################################################                                                                                                  
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=16"
set new3 = "  output_frequency=16"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=1.8"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 5, dt = 1.8...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5h_dt1p8.nc
echo 'Done computing tstep_type = 5, dt = 1.8...'

###########################################################                                                                                                        
##### tstep_type = 5, dt = 0.9, hydrostatic           #####                                                                                                        
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=32"
set new3 = "  output_frequency=32"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=0.9"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 5, dt = 0.9...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5h_dt0p9.nc
echo 'Done computing tstep_type = 5, dt = 0.9...'

###########################################################                                                                                                        
##### tstep_type = 5, dt = 0.45, hydrostatic          #####                                                                                                       
###########################################################                                                                                                        
set str1 = "tstep_type"
set str2 = "nmax"
set str3 = "output_frequency"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set linenum3 = `grep -n $str3 testnl.nl | cut -d: -f1`
set new1 = "  tsstep_type=5"
set new2 = "  nmax=64"
set new3 = "  output_frequency=64"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
sed -i "${linenum3}s,.*,$new3," testnl.nl
set str1 = "tsstep_type"
set str2 = "tstep"
set linenum1 = `grep -n $str1 testnl.nl | cut -d: -f1`
set linenum2 = `grep -n $str2 testnl.nl | cut -d: -f1`
set new1 = "  tstep_type=5"
set new2 = "  tstep=0.45"
sed -i "${linenum1}s,.*,$new1," testnl.nl
sed -i "${linenum2}s,.*,$new2," testnl.nl
echo 'Computing tstep_type = 5, dt = 0.45...'
mpirun -np $NCPU $EXEC < testnl.nl
mv $outputdir/dcmip2012_test31.nc $outputdir/tstep5h_dt0p45.nc
echo 'Done computing tstep_type = 5, dt = 0.45...'


module load anaconda3
cd $outputdir
cp $testdir/*.py $outputdir
python3 convergence_Th.py
python3 convergence_geo.py


echo message | ts
