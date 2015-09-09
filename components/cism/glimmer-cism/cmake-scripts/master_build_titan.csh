#!/bin/csh
# BETA Master build script for titan, last updated 01/17/2013 with v1699
# build the code in the 3* ways currently supported and submits some test jobs

# (1) *serial build not working yet, commented out
# (2) sometimes I had trouble at configure time with autoconf PGI build thinking it
# was in gnu mode, after having been in gnu mode. This was fixed by logging out then in again
# so something is not getting unloaded/swapped. Looking into it...

# (1) execute from the main seacism directory
# (2) set the next two commands 

#add logic at the top to decide which versions to build 
setenv TEST_DIR "/tmp/work/$USER/higher-order"
setenv CODE_DIR "/ccs/home/$USER/seacism"
cd $CODE_DIR
# 0 is a successful build
setenv build_no 0
# setting to 0 means don't build that version
setenv build_autoconf 1
setenv build_cmake 1

mkdir -p $TEST_DIR

# even if these are set in your env you need these when running the script
echo 'set the pgi env'
source /opt/modules/default/init/csh
module rm cmake netcdf python
module rm PrgEnv-pgi
module rm PrgEnv-gnu
module rm PrgEnv-intel
module rm PrgEnv-cray
module rm PrgEnv-pathscale
module add PrgEnv-pgi
module add netcdf-hdf5parallel/4.2.0 python subversion cmake
module swap xt-asyncpe xt-asyncpe/5.16
module list


# NEEDED AFTER A FRESH CHECKOUT
echo 'bootstrap'
./bootstrap >& bootstrap.out


###################################################################################
### AUTOCONF BUILDS
###################################################################################
echo "build_autoconf==$build_autoconf"
if ( $build_autoconf == 1 ) then

  # SERIAL BUILD WITH AUTOCONF PGI
  echo 'AUTOCONF PGI SERIAL BUILD'
  echo 'configure serial autoconf build'
  ./configure-scripts/titan-config-serial >& conf_autoconf_serial_pgi.out
  if ($status != 0) then
    echo 'ERROR: configure for autoconf, serial, PGI build'
    echo "cat $PWD/conf_autoconf_serial_pgi.out"
    @ build_no = 1
  endif
  echo 'make clean'
  make clean                              >& makeclean_autoconf_serial_pgi.out
  if ($status != 0) then
    echo 'ERROR: make clean for autoconf, serial, PGI build'
    echo "cat $PWD/makeclean_autoconf_serial_pgi.out"
    @ build_no = 1
  endif
  echo 'make serial'
  make                                    >& build_autoconf_serial_pgi.out 
  if ($status != 0) then
    echo 'ERROR: build for autoconf, serial, PGI build'
    echo "cat $PWD/build_autoconf_serial_pgi.out"
    @ build_no = 1
  endif

  if ( -e example-drivers/simple_glide/src/simple_glide ) then
    echo 'copy PGI autoconf serial executable to test directory'
    cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_serial
  endif
  # PARALLEL BUILD WITH AUTOCONF PGI
  echo 'AUTOCONF PGI PARALLEL BUILD'
  echo 'make distclean'
  make distclean                   >& makeclean_autoconf_parallel_pgi.out
  if ($status != 0) then
    echo 'ERROR: make clean for autoconf, parallel, PGI build'
    echo "cat $PWD/makeclean_autoconf_parallel_pgi.out"
    @ build_no = 1
  endif
  echo 'configure parallel autoconf build'
  ./configure-scripts/titan-config >& conf_autoconf_parallel_pgi.out
  if ($status != 0) then
    echo 'ERROR: configure for autoconf, parallel, PGI build'
    echo "cat $PWD/conf_autoconf_parallel_pgi.out"
    @ build_no = 1
  endif
  echo 'make parallel'
  make                             >& build_autoconf_parallel_pgi.out 
  if ($status != 0) then
    echo 'ERROR: build for autoconf, parallel, PGI build'
    echo "cat $PWD/build_autoconf_parallel_pgi.out"
    @ build_no = 1
  endif

  if ( -e example-drivers/simple_glide/src/simple_glide ) then
    echo 'copy PGI autoconf parallel executable to test directory'
    cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide
  endif

else # build with autoconf option
  echo 'not building autoconf code option'
endif # build with autoconf option



#Clean out anything left by autoconf builds
if ( $build_autoconf == 1 ) then
  cd $CODE_DIR
  make distclean >& distclean_postautoconf.out
endif



###################################################################################
### CMAKE BUILDS
###################################################################################
echo "build_cmake==$build_cmake"
if ( $build_cmake == 1 ) then

  # PARALLEL BUILD WITH CMAKE PGI
  echo 'CMAKE PGI PARALLEL BUILD'
  cd $CODE_DIR
  rm -rf titan-pgi
  mkdir titan-pgi
  cp cmake-scripts/titan-pgi-cmake titan-pgi
  cd $CODE_DIR/titan-pgi
  echo 'clean out the build dir'
  rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
  echo 'configure pgi cmake build'
  ./titan-pgi-cmake >& conf_cmake_parallel_pgi.out
  if ($status != 0) then
    echo 'ERROR: configure for cmake, parallel, PGI build'
    echo "cat $PWD/conf_cmake_parallel_pgi.out"
    @ build_no = 1
  endif
  echo 'make parallel pgi'
  make -j16         >& build_cmake_parallel_pgi.out
  if ($status != 0) then
    echo 'ERROR: build for cmake, parallel, PGI build'
    echo "cat $PWD/build_cmake_parallel_pgi.out"
    @ build_no = 1
  endif

  if ( -e example-drivers/simple_glide/src/simple_glide ) then
    echo 'copy pgi parallel executable to test directory'
    cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_pgi
  endif

  # PARALLEL BUILD WITH CMAKE GNU
  echo 'change to gnu env'
  module rm cmake netcdf-hdf5parallel/4.2.0 netcdf hdf5 python
  module swap PrgEnv-pgi PrgEnv-gnu
  module add cmake/2.8.6 python netcdf-hdf5parallel/4.2.0 
  module list
  echo 'CMAKE GNU PARALLEL BUILD'
  cd $CODE_DIR
  rm -rf titan-gnu
  mkdir titan-gnu
  cp cmake-scripts/titan-gnu-cmake titan-gnu
  cd $CODE_DIR/titan-gnu
  echo 'clean out the build dir'
  rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
  echo 'configure gnu cmake build'
  ./titan-gnu-cmake >& conf_cmake_parallel_gnu.out
  if ($status != 0) then
    echo 'ERROR: configure for cmake, parallel, GNU build'
    echo "cat $PWD/conf_cmake_parallel_gnu.out"
    @ build_no = 1
  endif
  echo 'make parallel gnu'
  make -j16         >& build_cmake_parallel_gnu.out
  if ($status != 0) then
    echo 'ERROR: build for cmake, parallel, GNU build'
    echo "cat $PWD/build_cmake_parallel_gnu.out"
    @ build_no = 1
  endif

  if ( -e example-drivers/simple_glide/src/simple_glide ) then
    echo 'copy gnu parallel executable to test directory'
    cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_gnu
  endif

else # build with cmake option
  echo 'not building cmake code option'
endif # build with cmake option

# execute tests on titan
# TODO the small jobs need to be combined into one ijob submission to get through the queue
if ($build_no == 1 ) then
  echo "no job sumbitted, build/builds failed"
else
  # simplest case, runs all builds and on a range of small processor counts 
  echo 'submitting jobs to compute nodes'
  #diagnostic dome test case
  #cd $TEST_DIR/reg_test/dome30/diagnostic
  #qsub ijob

  #evolving dome test case
  #cd $TEST_DIR/reg_test/dome30/evolving
  #qsub ijob

  # ISMIP test case A - not operational until BC set
  #cd $TEST_DIR/reg_test/ismip-hom-a/80km
  #qsub ijob

  # ISMIP test case C - not operational until BC set
  #cd $TEST_DIR/reg_test/ismip-hom-c/80km
  #qsub ijob

  # confined shelf to periodic BC
  #cd $TEST_DIR/reg_test/confined-shelf
  #qsub ijob

  # circular shelf to periodic BC
  #cd $TEST_DIR/reg_test/circular-shelf
  #qsub ijob

  # smaller GIS case to test realistic ice sheet configuration
  #cd $TEST_DIR/reg_test/gis_10km
  #qsub ijob

  # non regression test cases, default not run: 
  # large but not challenging case, to test large processor counts, not yet configured for titan
  #cd $TEST_DIR/dome500
  #qsub hopjob

  # high resolution GIS case to test realistic ice sheet configuration and longer time series, current setup gives
  #convergence problems
  #cd $TEST_DIR/gis_5km
  #qsub hopjob
endif

