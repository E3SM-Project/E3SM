## This is the automake file for building the Fortran examples
## for the PIO library.

# Ed Hartnett 7/17/19

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/set_flags.am

AM_FCFLAGS = -I$(top_srcdir)/src/flib

LDADD = ${top_builddir}/src/flib/libpiof.la	\
${top_builddir}/src/clib/libpioc.la

# Build the test for make check.
check_PROGRAMS = examplePio

examplePio_SOURCES = examplePio.f90

if RUN_TESTS
# Tests will run from a bash script.
TESTS = run_tests.sh
endif # RUN_TESTS

# Distribute the test script.
EXTRA_DIST = CMakeLists.txt run_tests.sh

# Clean up files produced during testing.
CLEANFILES = *.nc *.log *.mod
