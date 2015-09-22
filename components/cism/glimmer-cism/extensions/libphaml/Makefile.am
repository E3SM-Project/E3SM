FC = $(MPIFC)

if USE_PHAML_GRAPHICS
   PHAML_GRAPHICS = phaml_graphics
endif

bin_PROGRAMS = phaml_slave simple_phaml $(PHAML_GRAPHICS)

lib_LTLIBRARIES = libglide-phaml.la

libglide_phaml_la_LDFLAGS=-all-static
AM_LIBTOOLFLAGS="--tag=FC"
AM_FCFLAGS = $(NETCDF_FCFLAGS) $(FC_MODINC)$(top_builddir)/libglimmer-solve \
    $(FC_MODINC)$(top_builddir)/libglide $(FC_MODINC)$(top_builddir)/libglimmer \
    $(FC_MODINC)$(PHAML_PREFIX)/modules $(FC_MODINC)$(PHAML_PREFIX)/src

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak 
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)


libglide_phaml_la_SOURCES = phaml_user_mod.F90 phaml_support.F90 phaml_example.F90 \
                            phaml_pde.F90 phaml_example_pde.F90


#this slave program must be compiled for phaml to work
nodist_phaml_slave_SOURCES = $(PHAML_PREFIX)/src/slave.f90 
phaml_slave_LDADD = $(PHAML_PREFIX)/lib/libphaml.a \
    $(MPILIBS) libglide-phaml.la


#this builds the graphics slave
nodist_phaml_graphics_SOURCES = $(PHAML_PREFIX)/src/graphmain.f90
phaml_graphics_LDADD = $(PHAML_PREFIX)/lib/libphaml.a	\
        $(MPILIBS) libglide-phaml.la \
        -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lX11 -lXext -lm -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL

##############################################################
#a simple driver program
simple_phaml_SOURCES = simple_phaml.F90
simple_phaml_LDADD =  $(PHAML_PREFIX)/lib/libphaml.a \
    $(top_builddir)/libglide/libglide.la \
    $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
    $(top_builddir)/libglimmer/libglimmer-IO.la \
    $(top_builddir)/libglimmer/libglimmer.la \
    $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(MPILIBS) \
    libglide-phaml.la $(PHAML_PREFIX)/lib/libphaml.a 

f90_dep_sources = $(libglide_phaml_la_SOURCES)

-include f90_dep.mak
include $(top_srcdir)/extra_rules.am
