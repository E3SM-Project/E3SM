.. _unit-testing:

Fortran Unit Testing
====================

Introduction
------------

What is a unit test?
~~~~~~~~~~~~~~~~~~~~

A unit test is a fast, self-verying test of a small piece of code.
A single unit test typically covers 10s to 100s of lines of code (e.g., a single function or small module).
It typically runs in just milliseconds, and produces a simple pass/fail result.

Unit tests:

* Ensure that code remains correct as it is modified (in this respect, they complement the CIME system tests; the remaining bullet points are unique to unit tests).

* Ensure that new code is correct.

* Can help guide development, via test-driven development (TDD).

* Provide executable documentation of the intended behavior of a piece of code.

* Support development on your desktop machine.

Overview of unit test support in CIME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CIME comes with a set of tools to support building and running unit tests.
These consist of:

#. CMake tools to support building and running tests via CMake and CTest.

#. A python script that provides a simple front-end to the CMake-based tests.

The Fortran unit tests use `pFUnit <https://sourceforge.net/projects/pfunit/>`_, which is a Fortran unit testing framework that follows conventions of other xUnit frameworks (JUnit, etc.).

.. _running_unit_tests:

Running CIME's Fortran unit tests
---------------------------------

These instructions assume that you are using a machine that already has pFUnit installed, along with the necessary support in CIME.
If that is not the case, then see the section :ref:`adding_machine_support`.

From the top-level cime directory, you can run all of CIME's Fortran unit tests simply by running:

.. code-block:: shell

   > scripts/fortran_unit_testing/run_tests.py --build-dir MY_BUILD_DIR

where you can replace ``MY_BUILD_DIR`` with a path to the directory where you would like the unit test build files to be placed.
If you would like to ensure a completely clean build, you can use:

.. code-block:: shell

   > scripts/fortran_unit_testing/run_tests.py --build-dir `mktemp -d ./unit_tests.XXXXXXXX`

Once you have built the unit tests once (whether the build was successful or not), you can reuse the same build directory later to speed up the rebuild.
There are a number of useful arguments to ``run_tests.py``; for full usage information, run:

.. code-block:: shell

   > scripts/fortran_unit_testing/run_tests.py -h

If the build is successful, then you should get a message that looks like this:

.. code-block:: none

   ==================================================
   Running CTest tests for __command_line_test__/__command_line_test__.
   ==================================================

This will be followed by a list of tests, with a Pass/Fail message for each, like this:

.. code-block:: none

   Test project /Users/sacks/cime/unit_tests.0XHUkfqL/__command_line_test__/__command_line_test__
         Start  1: avect_wrapper
    1/17 Test  #1: avect_wrapper ....................   Passed    0.02 sec
         Start  2: seq_map
    2/17 Test  #2: seq_map ..........................   Passed    0.01 sec
         Start  3: glc_elevclass
    3/17 Test  #3: glc_elevclass ....................   Passed    0.01 sec   

You should then see a final message like this:

.. code-block:: none

   100% tests passed, 0 tests failed out of 17

For machines that have a serial build of pFUnit available for the default compiler, these unit tests are run automatically as part of ``scripts_regression_tests``.

.. _adding_machine_support:

How to add unit testing support on your machine
-----------------------------------------------

The below instructions assume that you have already ported CIME to your machine, by following the instructions in :doc:`/users_guide/porting-cime`.
Once you have done that, you can add unit testing support by building pFUnit on your machine and then pointing to the build in ``config_compilers.xml``.

At a minimum, you should do a serial build of pFUnit (without MPI or OpenMP), using the default compiler on your machine (according to ``config_machines.xml``).
That is the default used by ``run_tests.py``, and is required for ``scripts_regression_tests.py`` to run the unit tests on your machine.
Optionally, you can also provide pFUnit builds with other supported compilers on your machine.
If you'd like, you can also provide additional pFUnit builds with other combinations of MPI and OpenMP on or off.
However, at this time, no unit tests require parallel support, so there is no benefit gained by providing MPI-enabled builds.

Building pFUnit
~~~~~~~~~~~~~~~

To perform a serial build of pFUnit, follow these instructions:

#. Download pFUnit from https://sourceforge.net/projects/pfunit/

#. Set up your environment to be similar to the environment used in system builds of CIME.
   For example, load the appropriate compilers into your path.
   An easy way to achieve this is to run:

   .. code-block:: shell

      > $CIMEROOT/tools/configure --mpilib mpi-serial

   (with an optional ``--compiler`` argument; you'll also want to change the ``--mpilib`` argument if you're doing an MPI-enabled build).
   Then source either ``./.env_mach_specific.sh`` or ``./.env_mach_specific.csh``, depending on your shell.

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

If you'd like, you can then repeat this process with different compiler environments and/or different choices of ``-DMPI`` and ``-DOPENMP`` in the cmake step (each of these can have the value ``NO`` or ``YES``).
Make sure to choose a different installation directory for each of these, by setting the ``PFUNIT`` variable differently.

Adding to the xml file
~~~~~~~~~~~~~~~~~~~~~~

You then need to tell CIME about your pFUnit build(s).
To do this, specify the appropriate path(s) using the ``PFUNIT_PATH`` element in ``config_compilers.xml``.
For a serial build, this will look like:

.. code-block:: xml

     <PFUNIT_PATH MPILIB="mpi-serial" compile_threaded="false">$ENV{CESMDATAROOT}/tools/pFUnit/pFUnit3.2.8_cheyenne_Intel17.0.1_noMPI_noOpenMP</PFUNIT_PATH>

It is important that you provide the ``MPILIB`` and ``compile_threaded`` attributes.
``MPILIB`` should be ``mpi-serial`` for a pFUnit build with ``-DMPI=NO``, or the name of the mpi library you used for a pFUnit build with ``-DMPI=YES`` (e.g., ``mpich``; this should be one of this machine's MPI libraries specified by ``MPILIBS`` in ``config_machines.xml``).
``compile_threaded`` should be either ``true`` or ``false`` depending on the value of ``-DOPENMP``.

Once you have done this, you should be able to run the unit tests by following the instructions in :ref:`running_unit_tests`.

How to write a new unit test
----------------------------

TODO: Need to write this section.
This will draw on some of the information in sections 3 and 4 of https://github.com/NCAR/cesm_unit_test_tutorial (though without the clm and cam stuff).

General guidelines for writing unit tests
-----------------------------------------

Unit tests typically test a small piece of code (e.g., order 10 - 100 lines,
such as a single function or small-ish class).

Good unit tests are "FIRST"
(https://pragprog.com/magazines/2012-01/unit-tests-are-first):

* Fast (order milliseconds or less)

  * This means that, generally, they should not do any file i/o. Also, if you
    are testing a complex function, test it with a simple set of inputs - not a
    10,000-element array that will require a few seconds of runtime to process.

* Independent

  * This means that test Y shouldn't depend on some global variable that was
    created by test X. Dependencies like this cause problems if the tests run in
    a different order, if one test is dropped, etc.

* Repeatable

  * This means, for example, that you shouldn't generate random numbers in your
    tests.

* Self-verifying

  * This means that you shouldn't write a test that writes out its answers for
    manual comparison. Tests should generate an automatic pass/fail result.

* Timely

  * This means that the tests should be written *before* the production code
    (Test Driven Development), or immediately afterwards - not six months later
    when it's time to finally merge your changes onto the trunk, and have
    forgotten the details of what you have written. Much of the benefit of unit
    tests comes from developing them alongside the production code.

Good unit tests test a single, well-defined condition. This generally means that
you make a single call to the function / subroutine that you're testing, with a
single set of inputs. This means that you usually need multiple tests of the
function / subroutine, in order to test all of its possible behaviors. The main
reasons for testing a single condition in each test are:

* This makes it easier to pinpoint a problem when a test fails
* This makes it easier to read and understand the tests, allowing the tests to
  serve as useful documentation of how the code should operate

A good unit test has four distinct pieces:

#. **Setup**: e.g., create variables that will be needed for the routine you're
   testing. For simple tests, this piece may be empty.

#. **Exercise**: Call the routine you're testing

#. **Verify**: Call assertion methods to ensure that the results matched what
   you expected

#. **Teardown**: e.g., deallocate variables. For simple tests, this piece may be
   empty. **However, if this is needed, you should almost always do this
   teardown in the special tearDown routine, as discussed in the sections,**
   `Defining a test class in order to define setUp and tearDown methods`_ and
   `More on test teardown`_.

If you have many tests of the same subroutine, then you'll often find quite a
lot of duplication between the tests. It's good practice to extract major areas
of duplication to their own subroutines in the .pf file, which can be called by
your tests. This aids the understandability and maintainability of your
tests. pFUnit knows which subroutines are tests and which are "helper" routines
because of the ``@Test`` directives: You only add a ``@Test`` directive for your
tests, not for your helper routines.

More details on writing pFUnit-based unit tests
-----------------------------------------------

Assertion methods
~~~~~~~~~~~~~~~~~

pFUnit provides many assertion methods that you can use in the Verify step. Some
of the most useful are the following:

* ``@assertEqual(expected, actual)``

  * Ensures that expected == actual

  * Accepts an optional ``tolerance`` argument giving the tolerance for
    real-valued comparisons

* ``@assertLessThan(expected, actual)``

  * Ensures that expected < actual

* ``@assertGreaterThan(expected, actual)``

  * Ensures that expected > actual

* ``@assertLessThanOrEqual(expected, actual)``

* ``@assertGreaterThanOrEqual(expected, actual)``

* ``@assertTrue(condition)``

  * It's better to use the two-valued assertions above, if possible. For
    example, use ``@assertEqual(foo, bar)`` rather than ``@assertTrue(foo ==
    bar)``: the former gives more information if the test fails.

* ``@assertFalse(condition)``

* ``@assertIsFinite(value)``

  * Ensures that the result is not NaN or infinity

* ``@assertIsNan(value)``

  * Can be useful for failure checking, e.g., if your function returns NaN to
    signal an error

Comparison assertions accept an optional ``tolerance`` argument, which gives the
tolerance for real-valued comparisons.

In addition, all of the assertion methods accept an optional ``message``
argument, which gives a string that will be printed if the assertion fails. If
no message is provided, you will be pointed to the file and line number of the
failed assertion.

Defining a test class in order to define setUp and tearDown methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As noted in the comments in ``test_circle.pf``, the definition of a test class
(here, ``TestCircle``) is optional. However, it's convenient to define a minimal
test class when you first write a new .pf file:

.. code-block:: none

  @TestCase
  type, extends(TestCase) :: TestCircle
   contains
     procedure :: setUp
     procedure :: tearDown
  end type TestCircle

Defining this test class allows you to take advantage of some useful pFUnit
features like the setUp and tearDown methods.

If you define this test class, then you also need to:

* Define setUp and tearDown subroutines. These can start out empty:

  .. code-block:: Fortran

    subroutine setUp(this)
      class(TestCircle), intent(inout) :: this
    end subroutine setUp

    subroutine tearDown(this)
      class(TestCircle), intent(inout) :: this
    end subroutine tearDown

* Add an argument to each test subroutine, of class ``TestCircle`` (or whatever
  you called your test class). By convention, this argument is named ``this``.

Code in the setUp method will be executed before each test. This is convenient
if you need to do some setup that is the same for every test.

Code in the tearDown method will be executed after each test. This is often used
to deallocate memory. See the section, `More on test teardown`_ for details.

You can add any data or procedures to the test class. Adding data is
particularly useful, as this can be a way for the setUp and tearDown methods to
interact with your tests: The setUp method can fill a class variable with data,
which can then be used by your tests (accessed via
``this%somedata``). Conversely, if you want the tearDown method to deallocate a
variable, that variable cannot be local to your test subroutine. Instead, you
can make the variable a member of the class, so that the tearDown method can
access it.

So, for example, if you have this variable in your test class (as in the
example):

.. code-block:: Fortran

  real(r8), pointer :: somedata(:)

Then ``somedata`` can be created in the setUp method (if it needs to be the same
for every test). Alternatively, it can be created in each test routine that
needs it (if it differs from test to test, or some tests don't need it at
all). Its creation can look like:

.. code-block:: Fortran

  allocate(this%somedata(5))
  this%somedata(:) = [1,2,3,4,5]

Then your tearDown method can have code like this:

.. code-block:: Fortran

  if (associated(this%somedata)) then
    deallocate(this%somedata)
  end if

More on test teardown
~~~~~~~~~~~~~~~~~~~~~

All of the tests in a single test executable - which, for CIME, typically means
all of the tests defined in all ``.pf`` files in a single test directory - will
execute one after another in one run of the executable. This means that, if you
don't clean up after yourself, tests can interact with each other. In the best
case, this can mean you get a memory leak. In the worst case, it can mean that
the pass / fail status of tests depends on what other tests have run before
them, making your unit tests unrepeatable and unreliable. **As a general rule,
you should deallocate any pointers that your test allocated, reset any global
variables to some known, initial state, and do other, similar cleanup for
resources that may be shared by multiple tests.**

As described in the section, `Defining a test class in order to define setUp and
tearDown methods`_, code in the tearDown method will be executed after each
test. This is often used to do cleanup operations after each test. **Any
teardown like this should generally happen in this tearDown method. This is
because, if an assertion fails, the test aborts. So any teardown code in the
test method (following the failed assert statement) is skipped, which can lead
other tests to fail or give unexpected results. But this tearDown method is
still called in this case, making it a safe place to put teardown that needs to
be done regardless of whether the test passed or failed (which is the case for
most teardown).** In order for this to work, you sometimes need to move
variables that might otherwise be subroutine-local to the class - because the
tearDown method can access class instance variables, but not subroutine-local
variables.

Note that, in Fortran2003, allocatable variables are automatically deallocated
when they go out of scope, but pointers are not. So you need to explicitly
deallocate any pointers that have been allocated, either in test setup or in the
execution of the routine you're testing.

CIME makes extensive use of global variables: variables declared in some module,
which may be used (directly or indirectly) by the routine you're testing. If
your test has allocated or modified any global variables, it is important to
reset them to their initial state in the teardown portion of the
test. (Incidentally, this is just one of many reasons to prefer explicit
argument-passing over the use of global variables.)

pFUnit documentation and examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some pFUnit documentation is available here: http://pfunit.sourceforge.net/

If you download pFUnit (from
http://sourceforge.net/projects/pfunit/), you can find more
extensive documentation and examples in the following places. Among other
things, this can show you other assertion methods that are available:

* documentation/pFUnit3-ReferenceManual.pdf

* Examples/

* tests/

  * These are tests of the pFUnit code itself, written in pFUnit. You can see
    many uses of pFUnit features in these tests.


Finding more documentation and examples in CIME
-----------------------------------------------

Documentation of the unit test build system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CMake build infrastructure is in ``$CIMEROOT/src/externals/CMake``.

The infrastructure for building and running tests with ``run_tests.py`` is in
``$CIMEROOT/scripts/fortran_unit_testing``. That directory also contains some general
documentation about how to use the CIME unit test infrastructure (in the
``README`` file), and examples (in the ``Examples`` directory).

Finding more detailed examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At this point, there are many examples of unit tests in CIME, some simple and
some quite complex. You can find these by looking for files with the '.pf'
extension:

.. code-block:: shell

   > find . -name '*.pf'

You can also see examples of the unit test build scripts by viewing the
CMakeLists.txt files throughout the source tree.


