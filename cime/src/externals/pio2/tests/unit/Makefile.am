## This is the automake file for building the Fortran tests for the
## PIO library.

# Ed Hartnett 3/20/19

# Parallel builds don't currently work in this directory.
.NOTPARALLEL:

# Put together AM_CPPFLAGS and AM_LDFLAGS.
include $(top_srcdir)/set_flags.am

# Build the test for make check.
check_PROGRAMS = pio_unit_test_driver
pio_unit_test_driver_SOURCES = driver.F90
pio_unit_test_driver_LDADD = libglobal_vars.la libncdf_tests.la	\
libbasic_tests.la ${top_builddir}/src/flib/libpiof.la ${top_builddir}/src/clib/libpioc.la

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libglobal_vars.la libncdf_tests.la	\
libbasic_tests.la

# The convenience libraries depends on their source.
libglobal_vars_la_SOURCES = global_vars.F90
libncdf_tests_la_SOURCES = ncdf_tests.F90
libbasic_tests_la_SOURCES = basic_tests.F90

if RUN_TESTS
# Tests will run from a bash script.
TESTS = run_tests.sh
endif # RUN_TESTS

# Distribute the test script.
EXTRA_DIST = CMakeLists.txt run_tests.sh input.nl not_netcdf.ieee

# Clean up files produced during testing.
CLEANFILES = *.nc *.log *.mod
