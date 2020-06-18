03/2013 BFJ
02/2016 MT
10/2016 DMH
12/2019 MT  document some common cmake options

The CMAKE build system supports a number of user-configurable targets:
sweqx, preqx, preqx_acc, theta-l, swim, prim

Scripts which will CMake configure, build, construct namelists and run a simulation using these
targets, see:
  homme/test/sw_conservative/swtc[1256]ref.sh
  homme/test/jw_baroclinic/baro.job

A typical CMAKE command might look like:
cd $WDIR
cmake -C ~/e3sm/components/homme/cmake/machineFiles/edison.cmake \
    -DPREQX_PLEV=30 -DPREQX_NP=4 ~/acme/components/homme

The -C argument is machine specific settings. See the machineFile directory                   
for examples.  Most DOE machines are suported, as well as Darwin and RHEL.                    
                                                                                              
Some commonly used cmake options: (for both preqx and theta targets):                         
-DPREQX_USE_ENERGY=TRUE   add extra diagnostics to the log file                               
-DPREQX_USE_PIO=TRUE      turn on native grid output                                          
                          (default is interpolated lat/lon output)                            
                                                                                              
For performance testing, the model can be built without file output.  This                    
removes the need to link with PIO and netcdf.  See the astra.cmake machine                    
file and variable BUILD_HOMME_WITHOUT_PIOLIBRARY.                                             
                                                          
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

HOMME's regression tests use "cprnc" to compute differences between NETCDF files. HOMME
will attempt to build this from source (cime/tools/cprnc).  If you have trouble,
you can compile this independently and specify its location by setting
CPRNC_DIR in the cmake/machineFiles/yourmachine.cmake
To compile cprnc outside of CIME, edit the cprnc Makefile to comment out the Macros.cmake 
and set NETCDF_PATH, FC, and LDFLAGS=-L$(LIB_NETCDF) -lnetcdff -lnetcdf


DCMIP tests provide a standard means for testing and comparing the HOMME dycore with other dycores
both hydrostatic and nonhydrostatic. They have been placed in their own dcmip_test directory.
To run a DCMIP tests after configuring with cmake:
   1. navigate to the appropriate DCMIP working directory associated with the
   test case and model version, e.g. $WDIR/dcmip_tests/dcmip2016_test1_baroclinic_wave/theta-l
   2. "make install" to install test scripts and namelists.
   3. "./build.sh" to compile the executable
   4. edit one of the example jobscripts for your platform and then run it


The CMAKE code could use some cleanup. 
- user configured variables should not need to be prefixed by the exectuable name
  (i.e. -DNP=4, instead of -DPREQX_NP=4).  This will make the cmake code a lot simpler
  to maintain.
- the test cases are created after all the test executables are created. In keeping
  with cmake's tree-like directory approach, the tests should be associated with their
  test executable


