## This is a automake file, part of Unidata's netCDF package.
# Copyright 2005-2006, see the COPYRIGHT file for more information.

# Modified versions of nf_test routines used to test the use of
# F03 interface netcdf modules instead of netcdf include files

# This file builds and runs the following:

# The nf_test program, which tests the netCDF fortran 77 interface.

# The ftest program, which tests the netCDF fortran 77 V2
# interface. (Unless V2 is disabled).

# Test programs ftst_vars, ftst_vars2, ftst_vars3, ftst_vars4,
# ftst_vars5, ftst_types, ftst_types2, and ftst_groups, if netCDF-4 is
# enabled.


AM_FCFLAGS = -I$(top_srcdir)/fortran ${AM_CPPFLAGS}

# Some fortran compilers think your files should not have an .f90
# extension! The value of FCFLAGS_f90 is set in the configure script,
# based on the fortran compiler.
AM_FCFLAGS += $(FCFLAGS_f90) @MOD_FLAG@$(top_builddir)/fortran

# Need some stuff from ../nf_test
AM_FCFLAGS += -I$(top_srcdir)/nf_test

AM_FFLAGS = ${AM_FCFLAGS}

TESTS =

# All tests need to link to fortran and C libraries.
LDADD = ${top_builddir}/fortran/libnetcdff.la

# The create_fills shell creates the file fills.nc, needed by later
# tests.
TESTS += create_fills.sh
# Need access to top_srcdir
TESTS_ENVIRONMENT=TOPSRCDIR=${abs_top_srcdir}

# nf03_test is the main test program.
check_PROGRAMS = nf03_test
TESTS += nf03_test
nf03_test_SOURCES = module_tests.F90 f03lib_f_interfaces.f90 test03_get.F test03_put.F nf03_error.F nf03_test.F	\
test03_read.F test03_write.F util03.F f03lib.c

# programs that "use tests" from the tests.mod module depend on it being built first, needed for make -j
nf03_error.o nf03_test.o test03_get.o test03_put.o test03_read.o test03_write.o tst03_f77_v2.o util03.o: \
module_tests.o

# Cleaning up files created during the process.
CLEANFILES = scratch.nc test.nc copy.nc fills.nc

# Is the user building netCDF-4?
if USE_NETCDF4

# Add these netCDF-4 f77 test programs.
check_PROGRAMS += f03tst_vars f03tst_vars2 f03tst_vars3 f03tst_vars4 f03tst_vars5	\
f03tst_vars6 f03tst_types f03tst_types2 f03tst_types3 f03tst_groups
f03tst_types_SOURCES = module_tests.F90 f03tst_types.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_types2_SOURCES = module_tests.F90 f03tst_types2.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_types3_SOURCES = module_tests.F90 f03tst_types3.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_vars_SOURCES = module_tests.F90 f03tst_vars.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_vars2_SOURCES = module_tests.F90 f03tst_vars2.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_vars3_SOURCES = module_tests.F90 f03tst_vars3.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_vars4_SOURCES = module_tests.F90 f03tst_vars4.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_vars5_SOURCES = module_tests.F90 f03tst_vars5.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_vars6_SOURCES = module_tests.F90 f03tst_vars6.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
f03tst_groups_SOURCES = module_tests.F90 f03tst_groups.F f03lib_f_interfaces.f90 f03lib.c handle_err.F
TESTS += f03tst_vars f03tst_vars2 f03tst_vars3 f03tst_vars4 f03tst_vars5 \
         f03tst_vars6 f03tst_types f03tst_types2 f03tst_types3 f03tst_groups
CLEANFILES += ftst_vars.nc ftst_vars2.nc ftst_vars3.nc ftst_vars4.nc	\
              ftst_vars5.nc ftst_vars6.nc ftst_types.nc ftst_types2.nc \
              ftst_types3.nc ftst_groups.nc

# This is a netCDF-4 V2 test program.
if BUILD_V2
check_PROGRAMS += f03tst_v2
TESTS += f03tst_v2
CLEANFILES += ftst_v2.nc
f03tst_v2_SOURCES = f03tst_v2.F
endif # BUILD_V2

# Test parallel I/O
if TEST_PARALLEL
check_PROGRAMS += f03tst_parallel
CLEANFILES += ftst_parallel.nc
f03tst_parallel_SOURCES = f03tst_parallel.F
endif # TEST_PARALLEL
endif # USE_NETCDF4


# This is the fortran v2 test. It depends on the utilities being built
# to generate it's input file.
if BUILD_V2
check_PROGRAMS += f03test
TESTS += f03test
CLEANFILES += test.nc
f03test_SOURCES = module_tests.F90 f03test.F f03lib.c
endif # BUILD_V2

CLEANFILES += tests.mod

# test03_get.F and test03_put.f need to be distributed, so that the user
# need not have m4.
EXTRA_DIST = test03_get.F test03_put.F fills.cdl create_fills.sh	\
run_f77_par_test.sh CMakeLists.txt

test: check

# Did the user build the V2 F77 API? If so, run this test.
if BUILD_V2
check_PROGRAMS += tst03_f77_v2
TESTS += tst03_f77_v2
CLEANFILES += tst_f77_v2.nc temp.tmp
tst03_f77_v2_SOURCES = module_tests.F90 tst03_f77_v2.F
endif # BUILD_V2
