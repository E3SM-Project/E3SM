# BeTR regression test "meta-tests"

The test system needs to be the most bomb proof piece of code we
have. The 'mtest' directory contains the tests for the test system
itself.

## Running tests

If you make changes to the testing infracture, you *MUST* run the meta-test suite:

    # from the regression-tests directory, not mtest!
    make mtest
    # or
    python -m unittest discover --buffer --start-directory mtest/unit/


## Test Design

The testing infrastructure is designed to be heavily tested. Key principles are:

* read data into memory as soon as possible and isolate calles to the
  file system. Pass in-memory data objects to data processing routines
  to avoid working with files that are hard to test.

* no global variables

* no left justified code. all code must be in a subroutine so it can
  be isolated and tested, including main.

There are two types of tests for the test system, unit tests and
system tests. all are launched through a unit test framework.

The system tests are primarily expected failure tests. We seed the
test system with bad data and ensure that it catches it. system tests
should be run in 'check-only' mode when ever possible, or be extremely
short model runs.

## Code Coverage

A key aspect of being confident in the meta-test system is code
coverage. If you make changes to the test infrastructure, you must run
the code coverage tool and ensure that you are maintaining or
increasing code coverage before your changes will be accepted. Code
coverage reports are obtained by running the python coverage
tool. Install the coverage tool with:

    virtualenv env
    . env/bin/activate
    pip install coverage

Run the coverage report with:

    make test-coverage

from the regression-tests directory, not the mtest directory. This
will run the meta-tests in the coverage tool and output a simple text
summary to standard out, and create a detailed html coverage report:

    ...
    ...
    ...
    Name                                Stmts   Miss Branch BrPart  Cover
    ---------------------------------------------------------------------
    rtest_betr.py                         574    142    178     24    73%
    ---------------------------------------------------------------------
    TOTAL                                1127    142    178     24    84%
    
    5 files skipped due to complete coverage.
    For detailed analysis, see the html report at:
        coverage-html/index.html
    
