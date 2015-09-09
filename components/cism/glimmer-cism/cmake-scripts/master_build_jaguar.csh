#!/bin/csh
# Master build script for titan, last updated 1/3/2013 NOT YET TESTED
# build the code in the 3 ways currently supported and submits some test jobs
# note the serial build with autotools using SLAP is not supported on titan
# there are fewer tests run here than for jaguar since the allocation amount is smaller
 
# (1) execute from the main seacism directory
# (2) set the next two commands 

# user needs to set these two commands
setenv TEST_DIR "/tmp/work/4ue/SEACISM/cism_tests"
setenv CODE_DIR "$HOME/seacism"
cd $CODE_DIR
setenv build_no 0
setenv build_autoconf 1
setenv build_cmake 1
#mkdir -pv $TEST_DIR/configure_output

# even if these are set in your env you need these when running the script
echo 'set the pgi env'
module unload cmake netcdf python
module swap PrgEnv-gnu PrgEnv-pgi; module load cmake/2.8.7 python netcdf-hdf5-parallel/4.2.0 subversion 
#module load usg-default-modules/1.0

# NEEDED AFTER A FRESH CHECKOUT
echo 'bootstrap'
./bootstrap

echo $build_autoconf
if ($build_autoconf == 1 ) then
# PARALLEL BUILD WITH AUTOCONF PGI
echo 'make distclean'
make distclean
echo 'configure parallel autoconf build'
./configure-scripts/xk6-config >& conf_auto_pgi.out
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
echo 'flag to build cmake option:' $build_no
if ($build_cmake == 1 ) then
# set up directories for building cmake executables
rm -rf xk6-pgi
rm -rf xk6-gnu
mkdir xk6-pgi
mkdir xk6-gnu
cp cmake-scripts/jaguar-pgi-cmake xk6-pgi
cp cmake-scripts/jaguar-gnu-cmake xk6-gnu

# PARALLEL BUILD WITH CMAKE PGI
cd $CODE_DIR/xk6-pgi
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure pgi cmake build'
./jaguar-pgi-cmake >& conf_cmake_pgi.out
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
module swap PrgEnv-pgi PrgEnv-gnu; module load cmake/2.8.7 python netcdf-hdf5parallel/4.2.0
#module load usg-default-modules/1.0

cd $CODE_DIR/xk6-gnu
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure gnu cmake build'
./jaguar-gnu-cmake >& conf_gnu.out
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
# TODO the small jobs need to be combined into one ijob submission to get through the queue
if ($build_no == 1 ) then
  echo "no job sumbitted, build/builds failed"
else
# simplest case, runs all builds and on a range of small processor counts 
 echo 'submitting jobs to compute nodes'
#cd $TEST_DIR/dome30
#qsub ijob

# large but not challenging case, to test large processor counts, not yet configured for hopper
#cd $TEST_DIR/dome500
#qsub ijob

# ISMIP test case A - not activated until BC set
#cd $TEST_DIR/ismip-hom-a
#qsub ijob

# ISMIP test case C - not activated until BC set
#cd $TEST_DIR/ismip-hom-c
#qsub ijob

# confined shelf to periodic BC
#cd $TEST_DIR/confined-shelf
#qsub ijob

# circular shelf to periodic BC
#cd $TEST_DIR/circular-shelf
#qsub ijob

# smaller GIS case to test realistic ice sheet configuration
#cd $TEST_DIR/gis_10km
#qsub ijob

# high resolution GIS case to test realistic ice sheet configuration and longer time series, current setup gives
#convergence problems
#cd $TEST_DIR/gis_5km
#qsub ijob

endif
