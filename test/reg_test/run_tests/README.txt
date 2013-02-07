========================================================================
========================================================================
The scripts in this directory allow one to automate the building and running
running of tests.

========================================================================
========================================================================
Building:
========================================================================
The script build-all.sh manages the setup, configure and build of each test. 
Each test specifies information about the test in a .in file (eg. baro1a.in).
This information gets read producing a configuration script based on this 
information which builds the test executable out of place in the 
appropriate directory (HOMME_ROOT/build/preqx/baro1a/). The necessary 
files (namelist, reference solution, etc.) are also copied into the test 
executable directory. 

On Yellowstone, this script also creates lsf submission files for each test.

========================================================================
========================================================================
Running:
========================================================================
For now, the run_all.sh script only works on Yellowstone. Here, the lsf 
submission files produced by the build script are submitted to the queue and 
the results diffed against stdout results checked in to 
HOMME_ROOT/test/reg_test/results/yellowstone. These diff are concatenated 
into a single diff file for relatively simple viewing. 

========================================================================
========================================================================
Adding tests:
========================================================================
It is simple to add a test. Create a file your-test-name.in and define the
relevant information: test name, test flavor (eg. preqx, sweqx, etc.),
namelist, ncl files, vccord files etc. for your test. Then add the name of the
shell script "your-test-name.in" to the TESTS array in test-list.in. Be sure
to follow the existing examples so that the variable names are specified
correctly. 
Also, if you are adding a test to the repo, be certain to check that all of
the necessary files are "svn add"ed before you commit.

