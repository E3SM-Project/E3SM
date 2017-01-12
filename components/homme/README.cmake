03/2013 BFJ
02/2016 MT
10/2016 DMH

Please see the HOMME wiki for information on how to build HOMME using the CMake build system.
https://wiki.ucar.edu/display/homme/The+HOMME+CMake+build+and+testing+system

The CMAKE build system supports a number of user-configurable targets:
sweqx, preqx, preqx_acc, pese, swim, prim

Scripts which will CMake configure, build, construct namelists and run a simulation using these
targets, see:
  homme/test/sw_conservative/swtc[1256]ref.sh
  homme/test/jw_baroclinic/baro.job

A typical CMAKE command might look like:
cd $WDIR
cmake -C ~/acme/components/cmake/machineFiles/edison.cmake \
    -DPREQX_PLEV=30 -DPREQX_NP=4 ~/acme/components/homme

After running cmake with suitable command line options from a working directory WDIR,
it will create 

$WDIR/src/sweqx     directory containing user-configured sweqx executable
$WDIR/test_execs/swtcA     directory for swtcA executable
$WDIR/test_execs/swtcB     directory for swtcB executable
...and similarly for the preqx and preqx_acc targets.

$WDIR/utils         cprnc utility, and PIO and timing libraries
$WDIR/tests     directory containing all the HOMME regression tests

HOMME has a large regression test suite.  For instructions on running and adding
new tests, see homme/test/reg_test/README

DCMIP tests provide a standard means for testing and comparing the ACME HOMME dycore with other dycores
both hydrostatic and nonhydrostatic. They have been placed in their own dcmip_test directory for now.
To run a DCMIP tests, navigate to the appropriate directory and type make install to install test scripts and namelists.

The CMAKE code could use some cleanup. 
- user configured variables should not need to be prefixed by the exectuable name
  (i.e. -DNP=4, instead of -DPREQX_NP=4).  This will make the cmake code a lot simpler
  to maintain.
- the test cases are created after all the test executables are created. In keeping
  with cmake's tree-like directory approach, the tests should be associated with their


************************************************************************************************

***OBSOLETE***

03/2013 CGB and KJE and JER

HOMME now has a CMake build option for sweqx, swim, and preqx.
It is in the BETA test mode, alert Chris Baker or Kate Evans or Jennifer Ribbeckof problems/unclear instructions
(swim is the SW version of HOMME that uses trilinos in the implicit solve option)

1. mkdir BUILD_DIR somewhere, usually the main trunk directory
2. cp /bld/cmake-script/$APPROPRIATE_BUILD_SCRIPT into $BUILD_DIR 
3. export HOMME_ROOT="location_of trunk"
3. JAGUAR ONLY: export XTPE_LINK_TYPE='dynamic' needed right now TODO: put into script build process
3. Linux box ONLY: export Z_DIR='/usr/lib64' needed right now TODO: put into script build process
4. Modify script as appropriate (DEBUG or not etc)
 a. PLEV=# vertical levels
 b. NUM_POINTS=np
 c. -D CMAKE_INSTALL_PREFIX=where /bin/$EXE will sit 
5. ./$APPROPRIATE_BUILD_SCRIPT
5. make -j4
6. make install (where you told it to go)

most build failures are due to residual build info. In the build directory:
rm -rf CMakeCache.txt CMakeFiles src utils cmake_install.cmake Makefile bin install_manifest.txt *.h
to get a fresh build starting point.

