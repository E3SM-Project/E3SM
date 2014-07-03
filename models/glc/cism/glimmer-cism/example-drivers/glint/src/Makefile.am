bin_PROGRAMS  =  glint_example glex_ebm

include $(top_builddir)/libglimmer-trilinos/Makefile.Trilinos.export
TRILINOS_LIBS =  $(TRILINOS_LIBS_ALL)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = 
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

glint_example_SOURCES = glint_example.F90  trickLibtool.cpp
glex_ebm_SOURCES = glex_ebm.F90 glex_ebm_clim.F90   trickLibtool.cpp

LDADD = $(top_builddir)/libglint/libglint.la  $(top_builddir)/libglide/libglide.la \
        $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
        $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
        $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS) $(LAPACK_LIBS) $(PARDISO_LIBS) $(UMFPACK_LIBS) \
        $(TRILINOS_LIBS) $(FCLIBS) $(FLIBS)

AM_FCFLAGS = $(NETCDF_FCFLAGS) $(FC_MODINC)$(top_builddir)/libglimmer-solve $(FC_MODINC)$(top_builddir)/libglide \
             $(FC_MODINC)$(top_builddir)/libglimmer $(FC_MODINC)$(top_builddir)/libglint



f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -o $@ $^

# this fixes the build order: f90_dep.mak has .lo for eis_glide.o
# but since eis_glide isn't a library, there is no eis_glide.lo
# target, there is only a eis_glide.o target which never depends on
# the other targets, so it builds early and is missing it's .mod files
glex_ebm.o : glex_ebm.lo


-include f90_dep.mak
include $(top_srcdir)/extra_rules.am