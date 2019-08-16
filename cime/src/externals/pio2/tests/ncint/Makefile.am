## This is the automake file for building the netCDF integration layer
## tests.

# Ed Hartnett 7/3/19

# Put together AM_CPPFLAGS and AM_LDFLAGS.
AM_CPPFLAGS = -I$(top_srcdir)/src/clib
LDADD = ${top_builddir}/src/clib/libpioc.la

# Build the test for make check.
check_PROGRAMS = tst_pio_udf

if RUN_TESTS
# Tests will run from a bash script.
TESTS = run_tests.sh
endif # RUN_TESTS

# if RUN_TESTS
# # Tests will run from a bash script.
# TESTS = run_tests.sh
# endif # RUN_TESTS

# Distribute the test script.
EXTRA_DIST = run_tests.sh

# Clean up files produced during testing.
CLEANFILES = *.nc *.log
