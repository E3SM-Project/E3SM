
These directories contain working examples of the use of pFUnit to
perform unit testing.  They presume a working installation of pFUnit,
with the PFUNIT environment variable set to the directory where it is
installed.  

The GNUmakefile in this directory builds and runs the tests in serial
 or MPI mode. In general, to run the tests:
  - Set PFUNIT environment variable to directory where you installed pFUnit
and execute in the appropriate directory:
  - % make    # should build and run serial tests
  - % make MPI=YES   # should build and run all tests

Multiple approaches to unit-testing are possible with PFUNIT.

* Semi-automated (should always work) - developer indicates tests with
  semaphores in the source code and uses pFUnit preprocessing scripts
  to create a framework driver.  Assertions are done with macros that
  automatically track source file and line number, which improves
  clarity in several ways.  This approach reduces the effort to create
  and maintain tests but involves a bit of scripting "magic" that lies
  outside the domain of Fortran.  Real applications will likely need
  to adapt some aspects to their particular build system.

* Manual (Examples TBD) - developer manually inserts references to
  each test in the main driver.  Assertions are done with direct
  Fortran calls to pFUnit interfaces.  This approach is the most
  cumbersome, but also provides better understanding of how test
  assembly is managed.

* Automated (may never happen) - developer places all test code in a
  shared object library (.so).  pFUnit provides a static driver that
  can automatically execute tests in such a library given a regexp
  that indicates test naming conventions.  Default is "test*".


The current examples are:

* Simple - a very simple example that allows users to gain some
  confidence in the utility.  One test is purposefully failing to show what
  failure messages look like.

* Fixture - a more complicated example that involves subclassing TestCase.

* Halo - a very simple example showing how to create tests for MPI based
  procedures and execute them on varying processor counts.

* Parameterized - an example showing how to create and use a
  parameterized test case.

* MPI SimpleParameterized - a simple example showing how to create and use 
  a parameterized test case in MPI.

* Robust - demonstrates the use of the "-robust" option for the
  framework driver that allows for tests that crash or hang.

Note:
* These examples are not currently supported on Windows/CYGWIN, though a preliminary
  example is given in "Examples/Simple".
* For modifications and feature requests please see "sourceforge.net/projects/pfunit".


