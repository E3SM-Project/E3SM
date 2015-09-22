bin_PROGRAMS  =  eismint3_glide

include $(top_builddir)/libglimmer-trilinos/Makefile.Trilinos.export
TRILINOS_LIBS =  $(TRILINOS_LIBS_ALL)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak eismint3_vars.def
BUILT_SOURCES = eismint3_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

eismint3_glide_SOURCES = eismint3_glide.F90 eismint3_forcing.F90 eismint3_types.F90 eismint3_io.F90 trickLibtool.cpp
eismint3_glide_LDADD = $(top_builddir)/libglint/libglint.la  $(top_builddir)/libglide/libglide.la \
                       $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                       $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                       $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS) $(LAPACK_LIBS) $(PARDISO_LIBS) $(UMFPACK_LIBS)\
                       $(TRILINOS_LIBS) $(FCLIBS) $(FLIBS)

AM_FCFLAGS = $(NETCDF_FCFLAGS)  $(FC_MODINC)$(top_builddir)/libglimmer-solve $(FC_MODINC)$(top_builddir)/libglide \
             $(FC_MODINC)$(top_builddir)/libglimmer $(FC_MODINC)$(top_builddir)/libglint


#some special rules for automatically generating ncdf src code
eismint3_io.F90:	eismint3_vars.def $(top_srcdir)/libglimmer/ncdf_template.F90.in
			$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^

f90_dep.mak:	$(SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -o $@ $^

-include f90_dep.mak
# this fixes the build order: f90_dep.mak has .lo for eis_glide.o
# but since eis_glide isn't a library, there is no eis_glide.lo
# target, there is only a eis_glide.o target which never depends on
# the other targets, so it builds early and is missing it's .mod files
eismint3_glide.o : eismint3_glide.lo

include $(top_srcdir)/extra_rules.am
