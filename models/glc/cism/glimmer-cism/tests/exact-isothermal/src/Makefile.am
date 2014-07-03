bin_PROGRAMS = verif_glide
lib_LTLIBRARIES = libverif.la

include $(top_builddir)/libglimmer-trilinos/Makefile.Trilinos.export
TRILINOS_LIBS =  $(TRILINOS_LIBS_ALL)

variable_DATA = verif_vars.def
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
 
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = verif_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

libverif_la_SOURCES = verif.F90 verifBC.F90 verifD.F90 verif_io.F90
libverif_la_LDFLAGS=-all-static
verif_glide_SOURCES = verif_glide.F90 trickLibtool.cpp
verif_glide_LDADD = $(top_builddir)/libglide/libglide.la $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                    $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la libverif.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS) $(LAPACK_LIBS) $(PARDISO_LIBS) $(UMFPACK_LIBS) \
                    $(TRILINOS_LIBS) $(FCLIBS) $(FLIBS)

AM_FCFLAGS = $(NETCDF_FCFLAGS) $(FC_MODINC)$(top_builddir)/libglimmer-solve $(FC_MODINC)$(top_builddir)/libglide $(FC_MODINC)$(top_builddir)/libglimmer

#some special rules for automatically generating ncdf src code
verif_io.F90:	verif_vars.def $(top_srcdir)/libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^

-include f90_dep.mak

# this fixes the build order: f90_dep.mak has .lo for verif_glide.o
# but since verif_glide isn't a library, there is no verif_glide.lo
# target, there is only a verif_glide.o target which never depends on
# the other targets, so it builds early and is missing it's .mod files
verif_glide.o : verif_glide.lo

include $(top_srcdir)/extra_rules.am
