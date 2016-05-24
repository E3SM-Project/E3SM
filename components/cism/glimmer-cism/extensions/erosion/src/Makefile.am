lib_LTLIBRARIES = liberosion.la
check_PROGRAMS = test_rk4 test_geometry test_advect test_interpolate2d

variable_DATA = erosion_vars.def
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = erosion_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak *.mod

liberosion_la_SOURCES = erosion.F90 erosion_types.F90 erosion_setup.F90 erosion_io.F90 \
                       erosion_advect.F90 geometry.F90 runge_kutta.F90  erosion_transport.F90 \
                       glimmer_interpolate2d.F90 erosion_sediment.F90 advect_2ndmo.F90  erosion_transport_type.F90  

test_interpolate2d_SOURCES = test_interpolate2d.F90
test_geometry_SOURCES = test_geometry.F90
test_rk4_SOURCES = test_rk4.F90
test_advect_SOURCES = test_advect.F90

LDADD = liberosion.la $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                    $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)
AM_FCFLAGS = $(FC_MODINC)$(top_builddir)/libglide $(FC_MODINC)$(top_builddir)/libglimmer-solve $(FC_MODINC)$(top_builddir)/libglimmer \
             $(FC_MODINC)$(top_builddir)/tests/EISMINT/src $(FC_MODINC)$(top_builddir)/example-drivers/eis/src
test_advect_FCFLAGS = $(NETCDF_FCFLAGS)

#some special rules for automatically generating ncdf src code
erosion_io.F90:	erosion_vars.def $(top_srcdir)/libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^

-include f90_dep.mak
include $(top_srcdir)/extra_rules.am
