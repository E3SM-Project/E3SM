## This is the automake file for building the C examples for the PIO
## library.

# Ed Hartnett 5/7/18

# Link to our assembled library.
AM_LDFLAGS = ${top_builddir}/src/clib/libpio.la
AM_CPPFLAGS = -I$(top_srcdir)/src/clib

# Build the tests for make check.
check_PROGRAMS = example1 examplePio darray_no_async

# Tests will run from a bash script.
TESTS = run_tests.sh

# Distribute the test script.
EXTRA_DIST = run_tests.sh

# Clean up files produced during testing.
CLEANFILES = *.nc *.log
