.. _unit-testing:

Fortran Unit Testing
====================

Introduction
------------

What is a unit test?
~~~~~~~~~~~~~~~~~~~~

A unit test is a fast, self-verifying test of a small piece of code.
A single unit test typically covers 10s to 100s of lines of code; a single function or small module, for example.
It typically runs in milliseconds and produces a simple pass/fail result.

Unit tests:

* Ensure that code remains correct as it is modified. In this respect, unit tests complement the CIME system tests.

* Ensure that new code is correct.

* Can help guide development, via test-driven development (TDD).

* Provide executable documentation of the intended behavior of a piece of code.

* Support development on your desktop machine.

Overview of unit test support in CIME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CIME comes with a set of tools to support building and running unit tests.
These consist of:

#. CMake tools to support building and running tests via CMake and CTest.

#. A Python script that provides a simple front end for the CMake-based tests.

The Fortran unit tests use `pFUnit <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>`_, which is a Fortran testing framework that follows conventions of other xUnit frameworks.

.. _running_unit_tests:

Running CIME's Fortran unit tests
---------------------------------

These instructions assume that you are using a machine that already has pFUnit installed, along with the necessary support in CIME.
If that is not the case, see :ref:`adding_machine_support`.

From the top-level CIME directory, you can run all of CIME's Fortran unit tests by running:

.. code-block:: shell

   > scripts/fortran_unit_testing/run_tests.py --build-dir MY_BUILD_DIR

You can replace ``MY_BUILD_DIR`` with a path to the directory where you would like the unit test build files to be placed.
To ensure a completely clean build, use:

.. code-block:: shell

   > scripts/fortran_unit_testing/run_tests.py --build-dir `mktemp -d ./unit_tests.XXXXXXXX`

Once you have built the unit tests (whether the build was successful or not), you can reuse the same build directory later to speed up the rebuild.
There are a number of useful arguments to **run_tests.py**. For full usage information, run:

.. code-block:: shell

   > scripts/fortran_unit_testing/run_tests.py --help

If your build is successful, you will get a message like this:

.. code-block:: none

   ==================================================
   Running CTest tests for __command_line_test__/__command_line_test__.
   ==================================================

This will be followed by a list of tests, with a Pass/Fail message for each, like these examples:

.. code-block:: none

   Test project /Users/sacks/cime/unit_tests.0XHUkfqL/__command_line_test__/__command_line_test__
         Start  1: avect_wrapper
    1/17 Test  #1: avect_wrapper ....................   Passed    0.02 sec
         Start  2: seq_map
    2/17 Test  #2: seq_map ..........................   Passed    0.01 sec
         Start  3: glc_elevclass
    3/17 Test  #3: glc_elevclass ....................   Passed    0.01 sec

You will also see a final message like this:

.. code-block:: none

   100% tests passed, 0 tests failed out of 17

These unit tests are run automatically as part of **scripts_regression_tests** on machines that have a serial build of pFUnit available for the default compiler.

.. _adding_machine_support:

How to add unit testing support on your machine
-----------------------------------------------

The following instructions assume that you have ported CIME to your machine by following the instructions in :doc:`/users_guide/porting-cime`.
If you have done that, you can add unit testing support by building pFUnit on your machine and then pointing to the build in your **config_compilers.xml** file. Those processes are described in the following sections.

At a minimum, do a serial build of pFUnit (without MPI or OpenMP) using the default compiler on your machine.
That is the default that **run_tests.py** and that is required for **scripts_regression_tests.py** to run the unit tests on your machine.

Optionally, you can also provide pFUnit builds with other supported compilers on your machine.
You can also provide additional pFUnit builds with other combinations of MPI and OpenMP on or off.
At this time, however, no unit tests require parallel support so no benefit is gained by providing MPI-enabled builds.

Building pFUnit
~~~~~~~~~~~~~~~

For a serial build of pFUnit, follow these instructions:

#. Obtain pFUnit from https://github.com/Goddard-Fortran-Ecosystem/pFUnit (see
   https://github.com/Goddard-Fortran-Ecosystem/pFUnit#obtaining-pfunit for details; note
   that if you have an older version of cmake you may also need to use an older version of
   pFUnit)

#. Set up your environment to be similar to the environment used in CIME system builds.
   For example, load the appropriate compilers into your path.
   An easy way to achieve this is to run the following with an optional compiler argument:

   .. code-block:: shell

      > $CIMEROOT/tools/configure --mpilib mpi-serial

   If you are doing an MPI-enabled build, also change the ``--mpilib`` argument.
   Then source either **./.env_mach_specific.sh** or **./.env_mach_specific.csh**, depending on your shell.

   On some systems, you may still need to explicitly set the ``FC`` and ``CC`` environment
   variables, e.g., with:

   .. code-block:: shell

      > export FC=ifort
      > export CC=icc

#. For convenience, set the ``PFUNIT`` environment variable to point to the location where you want to install pFUnit. For example (in bash):

   .. code-block:: shell

      > export PFUNIT=/glade/p/cesmdata/cseg/tools/pFUnit/pFUnit3.2.8_cheyenne_Intel17.0.1_noMPI_noOpenMP

#. Configure and build pFUnit:

   .. code-block:: shell

      > mkdir build
      > cd build
      > cmake -DMPI=NO -DOPENMP=NO -DCMAKE_INSTALL_PREFIX=$PFUNIT ..
      > make -j 4

#. Run pFUnit's self-tests:

   .. code-block:: shell

      > make tests

#. Install pFUnit in the directory you specified earlier:

   .. code-block:: shell

      > make install

You can repeat this process with different compiler environments and/or different choices of ``-DMPI`` and ``-DOPENMP`` in the cmake step. (Each of them can have the value ``NO`` or ``YES``.)
Make sure to choose a different installation directory for each build by setting the ``PFUNIT`` variable differently.

Adding to the xml file
~~~~~~~~~~~~~~~~~~~~~~

After you build pFUnit, tell CIME about your build or builds.
To do this, specify the appropriate path(s) using the ``PFUNIT_PATH`` element in **config_compilers.xml**.
For a serial build, your setting will look like this example:

.. code-block:: xml

     <PFUNIT_PATH MPILIB="mpi-serial" compile_threaded="FALSE">$ENV{CESMDATAROOT}/tools/pFUnit/pFUnit3.2.8_cheyenne_Intel17.0.1_noMPI_noOpenMP</PFUNIT_PATH>

The ``MPILIB`` attribute should be either:

* ``mpi-serial`` for a pFUnit build where ``-DMPI=NO``, or

* the name of the MPI library you used for a pFUnit build where ``-DMPI=YES``. (For example, you might use ``mpich``, which should be one of the machine's MPI libraries specified by ``MPILIBS`` in **config_machines.xml**.)

The ``compile_threaded`` attribute should be either ``TRUE`` or ``FALSE`` depending on the value of ``-DOPENMP``.

Once you have specified the path for your build(s), you should be able to run the unit tests by following the instructions in :ref:`running_unit_tests`.

How to write a new unit test
----------------------------

.. todo:: Need to write this section. This will draw on some of the information in sections 3 and 4 of https://github.com/NCAR/cesm_unit_test_tutorial (though without the clm and cam stuff).

It should also introduce the role of .pf files, which are referenced several paragraphs later as if already explained.

General guidelines for writing unit tests
-----------------------------------------

Unit tests typically test a small piece of code, on the order of 10-100 lines, as in a single function or small class.

Good unit tests are **"FIRST"**:
(https://pragprog.com/magazines/2012-01/unit-tests-are-first):

* **Fast** (milliseconds or less). This means that, generally, they should not do any file i/o. Also, if you are testing a complex function, test it with a simple set of inputs rather than a 10,000-element array that will require a few seconds of runtime to process.

* **Independent**. This means that test Y shouldn't depend on some global variable that text X created. Such dependencies cause problems if the tests run in a different order, if one test is dropped, and so on.

* **Repeatable**. This means, for example, that you shouldn't generate random numbers in your tests.

* **Self-verifying**. Don't write a test that writes out its answers for manual comparison. Tests should generate an automatic pass/fail result.

* **Timely**. Write the tests *before* the production code (TDD) or immediately afterwards - not six months later when it's time to finally merge your changes onto the trunk and you have forgotten the details. Much of the benefit of unit tests comes from developing them concurrently with the production code.

Good unit tests test a single, well-defined condition. This generally means that
you make a single call to the function or subroutine that you're testing, with a
single set of inputs. Usually you need to run multiple tests in order to test
all of the unit's possible behaviors.

Testing a single condition in each test makes pinpointing problems easier when a test fails.
This also makes it easier to read and understand the tests, allowing them to serve as useful
documentation of how the code should operate.

A good unit test has four distinct pieces:

#. **Setup**: For example, creating variables that will be needed for the routine you're testing. For simple tests, this piece may be empty.

#. **Exercise**: Calling the routine you're testing.

#. **Verify**: Calling assertion methods (next section) to ensure that the results match what you expected.

#. **Teardown**: For example, deallocating variables. For simple tests, this piece may be empty. If it is needed, however, it is best done in the special tearDown routine discussed in `Defining a test class in order to define setUp and tearDown methods`_ and `More on test teardown`_.**

If you have many tests of the same subroutine, you may find quite a
lot of duplication. It's good practice to extract major areas of duplication to their own
subroutines in the **.pf** file, which your tests can call. This aids the understandability
and maintainability of your tests. pFUnit knows which subroutines are tests and which are
"helper" routines because of the ``@Test`` directives: You only add a ``@Test`` directive
for your tests, not for your helper routines.

More details on writing pFUnit-based unit tests
-----------------------------------------------

Assertion methods
~~~~~~~~~~~~~~~~~

pFUnit provides many assertion methods that you can use in the Verify step.
Here are some of the most useful:

=================================================    ===================================================================

``@assertEqual(expected, actual)``                   Ensures that expected == actual.
                                                     Accepts an optional ``tolerance`` argument giving the tolerance for
                                                     real-valued comparisons.

``@assertLessThan(expected, actual)``                Ensures that expected < actual.

``@assertGreaterThan(expected, actual)``             Ensures that expected > actual.

``@assertLessThanOrEqual(expected, actual)``

``@assertGreaterThanOrEqual(expected, actual)``

``@assertTrue(condition)``                           It is better to use the two-valued assertions above, if possible.
                                                     They provide more information if a test fails.

``@assertFalse(condition)``

``@assertIsFinite(value)``                           Ensures that the result is not NaN or infinity.

``@assertIsNan(value)``                              This can be useful for failure checking - for example, when your
                                                     function returns NaN to signal an error.

=================================================    ===================================================================

Comparison assertions accept an optional ``tolerance`` argument, which gives the
tolerance for real-valued comparisons.

All of the assertion methods also accept an optional ``message`` argument, which prints
a string if the assertion fails. If no message is provided, you will be pointed to the
file and line number of the failed assertion.

Defining a test class in order to define setUp and tearDown methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As noted in the comments in **test_circle.pf**, defining a test class is optional.
However, defining a minimal test class as shown here with ``TestCircle`` allows you
use some pFUnit features such as the setUp and tearDown methods.

.. code-block:: none

  @TestCase
  type, extends(TestCase) :: TestCircle
   contains
     procedure :: setUp
     procedure :: tearDown
  end type TestCircle

If you define this test class, you also need to:

* Define *setUp* and *tearDown* subroutines. These can start out empty:

  .. code-block:: Fortran

    subroutine setUp(this)
      class(TestCircle), intent(inout) :: this
    end subroutine setUp

    subroutine tearDown(this)
      class(TestCircle), intent(inout) :: this
    end subroutine tearDown

* Add an argument to each subroutine of the class. By convention, this argument is named ``this``.

Code in the setUp method is executed before each test. This is convenient
if you need to do some setup that is the same for every test.

Code in the tearDown method is executed after each test. This is often used
to deallocate memory. See `More on test teardown`_ for details.

You can add any data or procedures to the test class. Adding data is
particularly useful, as this can be a way for the setUp and tearDown methods to
interact with your tests: The setUp method can fill a class variable with data,
which your tests can then use (accessed via ``this%somedata``). Conversely, if
you want the tearDown method to deallocate a variable, the variable cannot be local
to your test subroutine. Instead, you make the variable a member of the class, so
that the tearDown method can access it.

Here is an example. Say you have this variable in your test class:

.. code-block:: Fortran

  real(r8), pointer :: somedata(:)

The setUp method can create ``somedata`` if it needs to be the same
for every test.

Alternatively, it can be created in each test routine that needs it if it
differs from test to test. (Some tests don't need it at all.) In that situation,
create it like this:

.. code-block:: Fortran

  allocate(this%somedata(5))
  this%somedata(:) = [1,2,3,4,5]

Your tearDown method then can have code like this:

.. code-block:: Fortran

  if (associated(this%somedata)) then
    deallocate(this%somedata)
  end if

More on test teardown
~~~~~~~~~~~~~~~~~~~~~

As stated in `Defining a test class in order to define setUp and tearDown methods`_,
code in the tearDown method is executed after each test, often to do cleanup operations.

Using the tearDown method is recommended because tests abort if an assertion fails.
The tearDown method is still called, however, so teardown that needs to be done
still gets done, regardless of pass/fail status. Teardown code might otherwise be
skipped, which can lead other tests to fail or give unexpected results.

All of the tests in a single test executable run one after another. For CIME, this
means all of the tests that are defined in all **.pf** files in a single test directory.

As a result, tests can interact with each other if you don't clean up after yourself.
In the best case, you might get a memory leak. In the worst case, the pass/fail status of tests
depends on which other tests have run previously, making your unit tests unrepeatable
and unreliable.

**To avoid this:**

* Deallocate any pointers that your test allocates.
* Reset any global variables to some known, initial state.
* Do other, similar cleanup for resources that are shared by multiple tests.

In Fortran2003, allocatable variables are deallocated automatically when they go
out of scope, but pointers are not. Explicitly deallocate any pointers that have
been allocated, either in test setup or in the execution of the routine
you are testing.

You might need to move some variables from subroutine-local to the class. This is
because the tearDown method can access class instance variables, but not subroutine-local
variables.

CIME makes extensive use of global variables that may be used directly or
indirectly by a routine you are testing. If your test has allocated or modified
any global variables, it is important to reset them to their initial state in the
teardown portion of the test.

Finding more documentation and examples
---------------------------------------

More detailed examples in CIME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are many examples of unit tests in CIME, some simple and some quite complex.
You can find them by looking for files with the ".pf" extension:

.. code-block:: shell

   > find . --name '*.pf'

You can also see examples of the unit test build scripts by viewing the
**CMakeLists.txt** files throughout the source tree.

Other pFUnit documentation sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extensive documentation and examples are included in the following when you obtain
pFUnit from https://github.com/Goddard-Fortran-Ecosystem/pFUnit:

* documentation/pFUnit3-ReferenceManual.pdf

* Examples/

* tests/

The tests are tests of the pFUnit code itself, written in pFUnit. They demonstrate
many uses of pFUnit features. Other documentation includes additional assertion
methods that are available.

Documentation of the unit test build system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CMake build infrastructure is in **$CIMEROOT/src/externals/CMake**.

The infrastructure for building and running tests with **run_tests.py** is in
**$CIMEROOT/scripts/fortran_unit_testing**. That directory also contains general
documentation about how to use the CIME unit test infrastructure (in the
**README** file) and examples (in the **Examples** directory).
