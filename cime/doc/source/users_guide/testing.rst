.. _testing:

**************************
Testing with create_test
**************************


**create_test** is the tool we use to test cime-driven models. It can be used as an easy way to
run a single basic test or an entire suite of tests.  create_test runs a test suite in parallel 
for improved performance.  create_test is the driver behind the automated nightly testing of 
cime-driven models.

Tests come in form::

  TESTTYPE.GRID.COMPSET[.MACHINE_COMPILER[.TESTMODS]]

Where TESTTYPE is how you want to do the test. Options for this value are listed in config_tests.xml.
A common test type would be ERS which stands for exact restart, which means you want to test the
grid and compset to see if they reproduce the same result upon restart as they would if no restart
took place. GRID and COMPSET are self-explanatory. MACHINE_COMPILER is optional; create_test will probe
the underlying machine and use the default compiler for that machine if this value is not
supplied. TESTMODS are optional and can be used to make some env XML changes prior to running
a test.

======================= =====================================================================================
TESTTYPE                Description
======================= =====================================================================================
   ERS                  Exact restart from startup (default 6 days + 5 days) 
                         Do an 11 day initial test - write a restart at day 6.    (file suffix: base) 

                         Do a 5 day restart test, starting from restart at day 6. (file suffix: rest) 

                         Compare component history files '.base' and '.rest' at day 11.
   ERS2                 Exact restart from startup  (default 6 days + 5 days).
                         Do an 11 day initial test without making restarts.     (file suffix: base)
 
                         Do an 11 day restart test stopping at day 6 with a restart, then resuming from restart at day 6. (file suffix: rest)
 
                         Compare component history files ".base" and ".rest" at day 11.
   ERT                  Exact restart from startup, default 2 month + 1 month (ERS with info DBUG = 1).
   IRT                  Exact restart from startup, (default 4 days + 7 days) with restart from interim file.
   ERIO                 Exact restart from startup with different PIO methods, (default 6 days + 5 days).
   ERR                  Exact restart from startup with resubmit, (default 4 days + 3 days).
   ERRI                 Exact restart from startup with resubmit, (default 4 days + 3 days). Tests incomplete logs option for st_archive.
   ERI                  hybrid/branch/exact restart test, default (by default STOP_N is 22 days) 

                        (1) ref1case

                            Do an initial run for 3 days writing restarts at day 3.

                            ref1case is a clone of the main case.

                            Short term archiving is on.

                        (2) ref2case (Suffix hybrid)

                            Do a hybrid run for default 19 days running with ref1 restarts from day 3,

                            and writing restarts at day 10. 

                            ref2case is a clone of the main case.

                            Short term archiving is on.

                        (3) case

                            Do a branch run, starting from restarts written in ref2case,

                            for 9 days and writing restarts at day 5.

                            Short term archiving is off.

                        (4) case (Suffix base)
       
                            Do a restart run from the branch run restarts for 4 days.

                            Compare component history files '.base' and '.hybrid' at day 19.

                            Short term archiving is off.
   ERP                  PES counts hybrid (OPENMP/MPI) restart bit for bit test from startup, (default 6 days + 5 days).
                         Initial PES set up out of the box

                         Do an 11 day initial test - write a restart at day 6.     (file suffix base)

                         Half the number of tasks and threads for each component.

                         Do a 5 day restart test starting from restart at day 6. (file suffix rest)

                         Compare component history files '.base' and '.rest' at day 11.

                         This is just like an ERS test but the tasks/threading counts are modified on restart
   PEA                  Single PE bit for bit test (default 5 days)
                         Do an initial run on 1 PE with mpi library.     (file suffix: base)

                         Do the same run on 1 PE with mpiserial library. (file suffix: mpiserial)

                         Compare base and mpiserial.
   PEM                  Modified PE counts for MPI(NTASKS) bit for bit test (default 5 days)
                         Do an initial run with default PE layout                               (file suffix: base)

                         Do another initial run with modified PE layout (NTASKS_XXX => NTASKS_XXX/2)  (file suffix: modpes)

                         Compare base and modpes
   PET                  Modified threading OPENMP bit for bit test (default 5 days)
                         Do an initial run where all components are threaded by default. (file suffix: base)

                         Do another initial run with NTHRDS=1 for all components.        (file suffix: single_thread)

                         Compare base and single_thread.

   PFS                  Performance test setup. (default 20 days)
   ICP                  CICE performance test.
   OCP                  POP performance test. (default 10 days)

   MCC                  Multi-driver validation vs single-instance. (default 5 days)
   NCK                  Multi-instance validation vs single instance - sequential PE for instances (default length)
                         Do an initial run test with NINST 1. (file suffix: base)

                         Do an initial run test with NINST 2. (file suffix: multiinst for both _0001 and _0002)

                         Compare base and _0001 and _0002.

   REP                  Reproducibility: Two identical runs are bit for bit. (default 5 days)
   SBN                  Smoke build-namelist test (just run preview_namelist and check_input_data).
   SMS                  Smoke startup test (default 5 days)
                         Do a 5 day initial test. (file suffix: base)
   SEQ                  Different sequencing bit for bit test. (default 10 days)
                         Do an initial run test with out-of-box PE-layout. (file suffix: base)

                         Do a second run where all root pes are at pe-0.   (file suffix: seq)

                         Compare base and seq.
   DAE                  Data assimilation test, default 1 day, two DA cycles, no data modification.
   PRE                  Pause-resume test: by default a bit for bit test of pause-resume cycling.
                         Default 5 hours, five pause/resume cycles, no data modification.

======================= =====================================================================================


Each test run by create_test includes the following mandatory steps:

* CREATE_NEWCASE: creating the create
* XML: xml changes to case based on test settings
* SETUP: setup case (case.setup)
* SHAREDLIB_BUILD: build sharedlibs
* MODEL_BUILD: build module (case.build)
* SUBMIT: submit test (case.submit)
* RUN: run test test

And the following optional phases:

* NLCOMP: Compare case namelists against baselines
* THROUGHPUT: Compare throughput against baseline throughput
* MEMCOMP: Compare memory usage against baseline memory usage
* MEMLEAK: Check for memleak
* COMPARE: Used to track test-specific comparions, for example, an ERS test would have a COMPARE_base_rest phase representing the check that the base result matched the restart result.
* GENERATE: Generate baseline results
* BASELINE: Compare results against baselines

Each test may be in one of the following states:

* PASS: The phase was executed successfully
* FAIL: We attempted to execute this phase, but it failed. If this phase is mandatory, no further progress will be made on this test. A detailed explanation of the failure should be in TestStatus.log.
* PEND: This phase will be run or is currently running but not complete

The current state of a test is represented in the file $CASEROOT/TestStatus

All output from the CIME infrastructure regarding this test will be put in the file $CASEROOT/TestStatus.log

A cs.status.$testid script will be put in the test root. This script will allow you to see the
current status of all your tests.

========================
Using create_test (E3SM)
========================

Usage will differ slightly depending on if you're using E3SM or CESM.

Using examples to illustrate common use cases

To run a test::

  ./create_test SMS.f19_f19.A

To run a test with a non-default compiler::

  ./create_test SMS.f19_f19.A --compiler intel

To run a test with baseline comparisons against baseline name 'master'::

  ./create_test SMS.f19_f19.A -c -b master

To run a test and update baselines with baseline name 'master'::

  ./create_test SMS.f19_f19.A -g -b master

To run a test with a non-default test-id::

  ./create_test SMS.f19_f19.A -t my_test_id

To run a test and use a non-default test-root for your case dir::

  ./create_test SMS.f19_f19.A -t $test_root

To run a test and use and put case, build, and run dirs all in the same root::

  ./create_test SMS.f19_f19.A --output-root $output_root

To run a test and force it to go into a certain batch queue::

  ./create_test SMS.f19_f19.A -q myqueue

To run a test and use a non-default project (can impact things like directory paths and acct for batch system)::

  ./create_test SMS.f19_f19.A -p myproj

To run two tests::

  ./create_test SMS.f19_f19.A SMS.f19_f19.B

To run a test suite::

  ./create_test e3sm_developer

To run a test suite excluding a specific test::

  ./create_test e3sm_developer ^SMS.f19_f19.A

See create_test -h for the full list of options

Interpreting test output is pretty easy, looking at an example::

  % ./create_test SMS.f19_f19.A

  Creating test directory /home/jgfouca/e3sm/scratch/SMS.f19_f19.A.melvin_gnu.20170504_163152_31aahy
  RUNNING TESTS:
    SMS.f19_f19.A.melvin_gnu
  Starting CREATE_NEWCASE for test SMS.f19_f19.A.melvin_gnu with 1 procs
  Finished CREATE_NEWCASE for test SMS.f19_f19.A.melvin_gnu in 4.170537 seconds (PASS)
  Starting XML for test SMS.f19_f19.A.melvin_gnu with 1 procs
  Finished XML for test SMS.f19_f19.A.melvin_gnu in 0.735993 seconds (PASS)
  Starting SETUP for test SMS.f19_f19.A.melvin_gnu with 1 procs
  Finished SETUP for test SMS.f19_f19.A.melvin_gnu in 11.544286 seconds (PASS)
  Starting SHAREDLIB_BUILD for test SMS.f19_f19.A.melvin_gnu with 1 procs
  Finished SHAREDLIB_BUILD for test SMS.f19_f19.A.melvin_gnu in 82.670667 seconds (PASS)
  Starting MODEL_BUILD for test SMS.f19_f19.A.melvin_gnu with 4 procs
  Finished MODEL_BUILD for test SMS.f19_f19.A.melvin_gnu in 18.613263 seconds (PASS)
  Starting RUN for test SMS.f19_f19.A.melvin_gnu with 64 procs
  Finished RUN for test SMS.f19_f19.A.melvin_gnu in 35.068546 seconds (PASS). [COMPLETED 1 of 1]
  At test-scheduler close, state is:
  PASS SMS.f19_f19.A.melvin_gnu RUN
    Case dir: /home/jgfouca/e3sm/scratch/SMS.f19_f19.A.melvin_gnu.20170504_163152_31aahy
  test-scheduler took 154.780044079 seconds

You can see that create_test informs the user of the case directory and of the progress and duration
of the various test phases.

===================
Managing baselines
===================

A big part of testing is managing your baselines (sometimes called gold results). We have provided
tools to help the user do this without having to repeat full runs of test cases with create_test.

bless_test_results: Takes a batch of cases of tests that have already been run and copy their
results to a baseline area.

compare_test_results: Takes a batch of cases of tests that have already been run and compare their
results to a baseline area.

Take a batch of results for the jenkins user for the testid 'mytest' and copy the results to
baselines for 'master'::

  ./bless_test_results -r /home/jenkins/e3sm/scratch/jenkins/ -t mytest -b master

Take a batch of results for the jenkins user for the testid 'mytest' and compare the results to
baselines for 'master'::

  ./compare_test_results -r /home/jenkins/e3sm/scratch/jenkins/ -t mytest -b master

=============
Adding tests
=============

Open the config/$model/tests.py file, you'll see a python dict at the top
of the file called _TESTS, find the test category you want to
change in this dict and add your testcase to the list.  Note the
comment at the top of this file indicating that you add a test with
this format: test>.<grid>.<compset>, and then there is a second
argument for mods.

========================
Scripts regression tests
========================

cime/scripts/tests/scripts_regression_tests.py is the suite of internal tests we run
for CIME. With no arguments, it will run the full suite. You can limit testing to a specific
test class or even a specific test within a test class.

Run full suite::

  ./scripts_regression_tests.py

Run a test class::

  ./scripts_regression_tests.py K_TestCimeCase

Run a specific test::

  ./scripts_regression_tests.py K_TestCimeCase.test_cime_case

If a test fails, the unittest module that drives scripts_regression_tests wil note the failure, but
won't print the output of the test until testing has completed. When there are failures for a
test, the case directories for that test will not be cleaned up so that the user can do a post-mortem
analysis. The user will be notified of the specific directories that will be left for them to
examine.
