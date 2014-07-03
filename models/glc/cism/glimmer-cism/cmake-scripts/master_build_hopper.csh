#!/bin/csh
# Master build script for hopper, last updated 11/09/2012 with v1617 
# build the code in the 4 ways currently supported and submits some test jobs
# there are fewer tests run here than for jaguar since the allocation amount is smaller
 
# (1) execute from the main seacism directory
# (2) set the next two commands 

#add logic at the top to decide which versions to build 
setenv TEST_DIR "$GSCRATCH/higher-order"
setenv CODE_DIR "$HOME/seacism"
cd $CODE_DIR
setenv build_no 0
setenv build_autoconf 1
setenv build_cmake 1
#mkdir -pv $TEST_DIR/configure_output

# even if these are set in your env you need these when running the script
echo 'set the pgi env'
module unload cmake netcdf python
module swap PrgEnv-gnu PrgEnv-pgi; module load cmake/2.8.7 python netcdf-hdf5parallel/4.2.0 subversion usg-default-modules/1.0

# NEEDED AFTER A FRESH CHECKOUT
echo 'bootstrap'
./bootstrap

echo $build_autoconf
if ($build_autoconf == 1 ) then
#SERIAL BUILD WITH AUTOCONF
echo 'make distclean'
make distclean 
echo 'configure serial autoconf build'
./configure-scripts/hopper-config-serial >& conf_serial.out
echo 'make serial'
make >& serial_build.out 

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy serial executable to test directory'
 cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_serial
else
 echo "autoconf parallel build failed, no executable"
 @ build_no = 1
endif

# PARALLEL BUILD WITH AUTOCONF PGI
echo 'make distclean'
make distclean
echo 'configure parallel autoconf build'
./configure-scripts/hopper-config-cesmtimers >& conf_auto_pgi.out
echo 'make parallel'
make >& auto_pgi_build.out 

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy autoconf parallel executable to test directory'
 cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide
else
 echo "autoconf parallel build failed, no executable"
 @ build_no = 1
endif

else # build with autoconf option
 echo 'not building autoconf code option'
endif # build with autoconf option

echo $build_no
echo 'flag to build cmake option:' $build_cmake
if ($build_cmake == 1 ) then
# set up directories for building cmake executables
rm -rf xe6-pgi
rm -rf xe6-gnu
mkdir xe6-pgi
mkdir xe6-gnu
cp cmake-scripts/hopper-pgi-cmake-cesmtimers xe6-pgi
cp cmake-scripts/hopper-gnu-cmake-cesmtimers xe6-gnu

# PARALLEL BUILD WITH CMAKE PGI
cd $CODE_DIR/xe6-pgi
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure pgi cmake build'
./hopper-pgi-cmake-cesmtimers >& conf_cmake_pgi.out
echo 'make parallel pgi'
make -j 4 >& cmake_pgi_build.out

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy pgi parallel executable to test directory'
 cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_pgi
else
 echo "cmake pgi build failed, no executable"
 @ build_no = 1
endif

cd $CODE_DIR
make distclean

echo $build_no
# PARALLEL BUILD WITH CMAKE GNU
echo 'change to gnu env'
module unload cmake netcdf-hdf5parallel/4.2.0 python
module swap PrgEnv-pgi PrgEnv-gnu; module load cmake/2.8.7 python netcdf-hdf5parallel/4.2.0 usg-default-modules/1.0

cd $CODE_DIR/xe6-gnu
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure gnu cmake build'
./hopper-gnu-cmake-cesmtimers >& conf_gnu.out
echo 'make parallel gnu'
make -j 4 >& cmake_gnu_build.out

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy gnu parallel executable to test directory'
 cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_gnu
else
 echo "cmake gnu build failed, no executable"
 @ build_no = 1
endif

else # build with cmake option
 echo 'not building cmake code option'
endif # build with cmake option

echo $build_no
# execute tests on hopper
# TODO the small jobs need to be combined into one hopjob submission to get through the queue
if ($build_no == 1 ) then
  echo "no job sumbitted, build/builds failed"
else
# simplest case, runs all builds and on a range of small processor counts 
 echo 'submitting jobs to compute nodes'
#diagnostic dome test case
cd $TEST_DIR/reg_test/dome30/diagnostic
qsub hopjob

#evolving dome test case
cd $TEST_DIR/reg_test/dome30/evolving
qsub hopjob

# ISMIP test case A - not operational until BC set
cd $TEST_DIR/reg_test/ismip-hom-a/80km
qsub hopjob

# ISMIP test case C - not operational until BC set
cd $TEST_DIR/reg_test/ismip-hom-c/80km
qsub hopjob

# confined shelf to periodic BC
cd $TEST_DIR/reg_test/confined-shelf
qsub hopjob

# circular shelf to periodic BC
cd $TEST_DIR/reg_test/circular-shelf
qsub hopjob

# smaller GIS case to test realistic ice sheet configuration
cd $TEST_DIR/reg_test/gis_10km
qsub hopjob

# non regression test cases, default not run: 

# large but not challenging case, to test large processor counts, not yet configured for hopper
#cd $TEST_DIR/dome500
#qsub hopjob

# high resolution GIS case to test realistic ice sheet configuration and longer time series, current setup gives
#convergence problems
#cd $TEST_DIR/gis_5km
#qsub hopjob

endif
