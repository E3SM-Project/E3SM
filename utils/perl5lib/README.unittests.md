Unit test directory for CIME perl5lib
=====================================

The goal of the unit test suite for cime perl5lib is to have machine
independent tests that will run anywhere in seconds and have high
coverage of the critical decision making code.

Running
-------

To run the unit tests for cime perl5lib:

    $ prove test_p5lib.pl 

The test output has been optimized to most useful when run through
prove (a standard part of all modern perl distributions). Developers
who prefer standard Test::More output can run the tests by running the
script:

    $ ./test_p5lib.pl


Creating New Tests
------------------

To create a new test suite, copy the file:

t/template_test_XXX.pm 

to t/test_what_you_want_to_test.pm

Note: the contents of template_test_XXX.pm are valid tests and are
part of the test suite. This file should always be valid.

New test files are picked up automatically by test_p5lib.pl.

Inside the new test module:

* WWW, XXX, YYY, ZZZ are used as place holder text for things that
  need to be replaced.

* startup and shutdown - common fixtures for all tests. These methods
  are only called once for each suite. The objects in these functions
  should NOT BE MODIFIED by any tests, e.g. config_cache,
  namelist_definition.xml, namelist_defaults.xml

* setup and teardown - common fixtures for all tests. These methods
  are called once for each test. Objects that ARE MODIFIED should be
  created here. For example, the namelist object should be created
  here so each test starts with a clean object.

* Tests are automatically detected by functions starting with
  'test_'. The proposed naming conventions are 'test_XXX__YYY' where
  XXX is the function or namelist variable being tested. YYY is a
  descriptor of the branch (if phys==clm4_5) or condition being
  tested. This will make the failure output more useful.

* To create each test, simply create a minimal set of input opts and
  pass them into the function you want to test, then check the
  results.

* Multiple tests can be placed in each function (including loops), but
  this make it harder to understand what failed and why. By using one
  test per test function, the failure messages will be more useful.

* Rather than describing the test in comments, it should be described
  in the message string. The message string should be printed to the
  screen when the test fails with "... || diag($msg)". This provides a
  useful failure message without duplicating information in comments
  that are probably out of date.

* To keep the tests machine independent, we don't want to require the
  presence of the model input data repository or other files at system
  specific locations. There are two options meet this requirement but
  still allow for testing of key functionality:

  1. Write and destroy small temporary input files as part of the test
  fixtures. This keeps the tests self contained and focused. Any test
  that simply tests for the existance of a file should not include a
  large source file in the repo. Just create an empty text
  file. (**prefered**)

  2. Modify or add new mock xml files in cime/utils/perl5lib/t/input
  and point to them as part of the tets. (**discouraged**)

  Side effects of this convention is that we minimize hard code paths
  and ensure we aren't depending on any one projects conventions.

