lib_LTLIBRARIES = libglint.la

variable_DATA = glint_vars.def glint_mbal_vars.def
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
 
EXTRA_DIST = f90_dep.mak 
BUILT_SOURCES = glint_io.F90 glint_mbal_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

libglint_la_SOURCES = glint_daily_pdd.F90 glint_pdd.F90 glint_constants.F90 \
                      glint_type.F90 glint_initialise.F90 glint_io.F90 \
                      glint_timestep.F90 glint_main.F90 \
                      glint_interp.F90 glint_global_interp.F90 \
                      glint_global_grid.F90 glint_routing.F90 \
                      glint_precip_param.F90 glint_ebm.F90 \ 
                      glint_mbal.F90 glint_mbal_coupling.F90 glint_mbal_io.F90 \
                      glint_mpinterp.F90 glint_commandline.F90 \
                      glint_example_clim.F90
libglint_la_LDFLAGS=-all-static

AM_FCFLAGS = -I$(top_srcdir)/libglide $(FC_MODINC)../libglimmer $(FC_MODINC)../libglide $(FC_MODINC)../libglimmer-solve $(NETCDF_FCFLAGS)

#some special rules for automatically generating ncdf src code
glint_io.F90:	glint_vars.def ../libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^
glint_mbal_io.F90:	glint_mbal_vars.def ../libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^


-include f90_dep.mak

glide_mask.inc : $(top_srcdir)/libglide/glide_mask.inc

include $(top_srcdir)/extra_rules.am
