#!/bin/csh



# Master build script for mac laptops. Last updated 2/28/2013 by SFP.
# This is a hacked version of Kate's original script for use on Hopper.
# For now, only supports parallel build with Trilinos using gnu and cmake. 
# Only a subset of the small, standard tests are run, on both 1 and 4 procs.
 
# (1) execute from the builds/blizzard-gnu subdirectory of CISM

#add logic at the top to decide which versions to build 

# PARALLEL BUILD WITH CMAKE 

# setenv TEST_DIR "/USERS/$USER/work/modeling/cism/seacism-oceans11/tests/higher-order"

# 5/7/2014 DMR -- added performance tests:

## This will automatically submit dome60-500 ijobs. gis_1km and gis_4km will not be submitted
## automatically because you will have to build and run Felix/Albany on hopper first. Once you do that,
## you can go to lines #193-194, 197-198, 201-202, and uncomment them.
setenv PERF_TEST 0

@ run_perf_tests = (($1 == run-perf-tests) || ($2 == run-perf-tests) || ($3 == run-perf-tests) || ($4 == run-perf-tests) || ($5 == run-perf-tests))

if ($run_perf_tests) then
  setenv PERF_TEST 1
endif

@ skip_build_set = (($1 == skip-build) || ($2 == skip-build) || ($3 == skip-build) || ($4 == skip-build) || ($5 == skip-build))

@ no_copy_set = (($1 == no-copy) || ($2 == no-copy) || ($3 == no-copy) || ($4 == no-copy) || ($5 == no-copy))

@ skip_tests_set = (($1 == skip-tests) || ($2 == skip-tests) || ($3 == skip-tests) || ($4 == skip-tests) || ($5 == skip-tests))

#**!move this and source it to your .bashrc (wherever your higher-order directory is located)
#setenv TEST_DIR /lustre/atlas/scratch/$USER/cli062/higher-order

if (! -d $TEST_DIR) mkdir -p $TEST_DIR

setenv TEST_SUITE_DEFAULT_LOC  http://oceans11.lanl.gov/cism/livv
#setenv TEST_SUITE_DEFAULT_LOC /ccs/proj/cli062/test_suite

setenv build_problem 0

set COMPILER_NAME = gnu
set PLATFORM_NAME = blizzard

# set PLATFORM_NAME = $1
# set COMPILER_NAME = $2

set CMAKE_SCRIPT = $PLATFORM_NAME'-'$COMPILER_NAME'-cmake'
set CMAKE_CONF_OUT = 'conf_'$COMPILER_NAME'.out'
set CMAKE_BUILD_OUT = 'cmake_'$COMPILER_NAME'_build.out'
#set CISM_RUN_SCRIPT = $PLATFORM_NAME'job' 
#set CISM_RUN_SCRIPT = 'hopjob'
set CISM_RUN_SCRIPT = 'ijob_linux' 
set CISM_VV_SCRIPT = $PLATFORM_NAME'_VV.bash'
#set CISM_VV_SCRIPT = 'rhea_VV.bash'

echo
echo 'To use this script, type: csh '$PLATFORM_NAME'-'$COMPILER_NAME'-build-and-test.csh'
echo
#echo 'For a quick test (dome only), type: csh '$PLATFORM_NAME'-'$COMPILER_NAME'-build-and-test.csh quick-test'
echo
echo "Call with no-copy to prevent copying of the reg_test and livv defaults."
echo "Call with run-perf-tests to run the performance tests."
echo "Call with skip-tests to skip testing (builds executable and copies it to TEST_DIR)."


echo
echo 'See the LIVV documentation for instructions on setting up the test directory (TEST_DIR).'
echo


#echo 'The following environment variables must be set: TEST_DIR, GLIMMER_TRILINOS_DIR'
#echo 'Examples (place in .cshrc or .bashrc):'
#echo 'csh, tcsh:  setenv GLIMMER_TRILINOS_DIR "/Users/$USER/Trilinos/gcc-build/install"'
#echo 'bash:       export GLIMMER_TRILINOS_DIR="/Users/$USER/Trilinos/gcc-build/install"'
echo
echo 'Setting TEST_DIR to the location: '
echo 'TEST_DIR =' $TEST_DIR
echo 'TEST_DIR must also be set in your .bashrc file.'

# PARALLEL BUILD WITH CMAKE


if ($skip_build_set == 0) then

echo
echo "Configuring and building in directory: " $PWD
echo 

echo 'Configuring '$COMPILER_NAME' cmake build...'
source ./$CMAKE_SCRIPT >& $CMAKE_CONF_OUT
echo 'Making parallel '$COMPILER_NAME'...'
make -j 8 >& $CMAKE_BUILD_OUT

#if ( -e example-drivers/simple_glide/src/simple_glide ) then
# echo 'Copying '$COMPILER_NAME' parallel simple_glide_'$COMPILER_NAME' to test directory'
# cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_$COMPILER_NAME
#else
# echo "cmake '$COMPILER_NAME' build failed, no executable"
# @ build_problem = 1
#endif

if ( -e cism_driver/cism_driver ) then
 echo 'Copying '$COMPILER_NAME' parallel cism_driver_'$COMPILER_NAME' to test directory'
 cp -f cism_driver/cism_driver $TEST_DIR/cism_driver_$COMPILER_NAME
else
 echo "cmake '$COMPILER_NAME' build failed, no executable"
 @ build_problem = 1
endif

endif # skip_build_set

if ($build_problem == 1) then
  echo "No job submitted -- cmake build failed."
else  # execute tests:
 
 # Make copy of test suite in $TEST_DIR:
if (! ($no_copy_set)) then
 echo "Copying default reg_test and LIVV to $TEST_DIR"
 pushd . > /dev/null
 cd $TEST_DIR
 if ( -e reg_test_default.tgz ) rm -f reg_test_default.tgz 
 wget $TEST_SUITE_DEFAULT_LOC/reg_test_default.tgz
 tar xfz reg_test_default.tgz
 popd > /dev/null

 if ($PERF_TEST) then
    echo "Copying default perf_test to $TEST_DIR"
   pushd . > /dev/null
   cd $TEST_DIR
   if ( -e perf_test_default.tgz ) rm -f perf_test_default.tgz 
   wget $TEST_SUITE_DEFAULT_LOC/perf_test_default.tgz
   tar xfz perf_test_default.tgz
   popd > /dev/null
 endif

 cp -rf ../../tests/higher-order/livv $TEST_DIR
endif

if ($skip_tests_set) then
   echo "Skipping tests."
   exit
endif

csh $TEST_DIR/livv/run_livv_default_tests.csh $TEST_DIR $CISM_RUN_SCRIPT $PERF_TEST $CISM_VV_SCRIPT
echo "Back in build-and-test script, exiting."
exit


