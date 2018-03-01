## This is an automake file, part of Unidata's netCDF package.
# Copyright 2005-2011, see the COPYRIGHT file for more information.

# This file controls the building of the c library
# that interfaces between fortran and C.

# The (temporary) library we are building.
# The version number is specified using the following
# guideline:
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
noinst_LTLIBRARIES = libnetcdfc.la
libnetcdfc_la_LDFLAGS = $(AM_LDFLAGS) -version-number 6:0:1

# Tell make to install these as a header file.
include_HEADERS = nfconfig.inc

# Turn on some special flags for building DLL.
AM_CPPFLAGS = -DDLL_EXPORT
if BUILD_DLL
AM_CPPFLAGS += -DDLL_EXPORT
libnetcdfc_la_LDFLAGS += -no-undefined -Wl,--output-def,netcdffdll.def
endif # BUILD_DLL

# These are the source files for the fortran library.
libnetcdfc_la_SOURCES = fort-attio.c fort-control.c fort-dim.c		\
fort-genatt.c fort-geninq.c fort-genvar.c fort-lib.c fort-misc.c	\
fort-v2compat.c fort-vario.c fort-var1io.c fort-varaio.c fort-varmio.c	\
fort-varsio.c cfortran.h fort-lib.h ncfortran.h

# Is the user building netCDF-4?
if USE_NETCDF4
# All netCDF-4 F77 functions are in one code file.
libnetcdfc_la_SOURCES += fort-nc4.c
endif

## The include file and man page are both built by make.
#BUILT_SOURCES = netcdf.inc
#
#netcdf.inc: netcdf2.inc netcdf3.inc netcdf4.inc nfconfig.inc
#	echo '!     NetCDF-3.' > netcdf.inc
#	cat @top_srcdir@/fortran/netcdf3.inc >> netcdf.inc
#if BUILD_V2
#	echo >> netcdf.inc
#	echo '!     NetCDF-2.' >> netcdf.inc
#	cat @top_srcdir@/fortran/netcdf2.inc >> netcdf.inc
#endif
#if USE_NETCDF4
#	echo >> netcdf.inc
#	echo '!     NetCDF-4.' >> netcdf.inc
#	cat @top_srcdir@/fortran/netcdf4.inc >> netcdf.inc
#endif
#if USE_LOGGING
#	echo >> netcdf.inc
#	echo '!     This is to turn on netCDF internal logging.' >> netcdf.inc
#	echo '      integer nf_set_log_level' >> netcdf.inc
#	echo '      external nf_set_log_level' >> netcdf.inc
#endif

# The F77 API man page.
man_MANS = netcdf_f77.3

# This rule generates the F77 manpage.
ARGS_MANPAGE = -DAPI=FORTRAN
if USE_NETCDF4
ARGS_MANPAGE += -DNETCDF4=TRUE
endif
if BUILD_DAP
ARGS_MANPAGE += -DDAP=TRUE
endif
if BUILD_PARALLEL
ARGS_MANPAGE += -DPARALLEL_IO=TRUE
endif
netcdf_f77.3: $(top_srcdir)/docs/netcdf.m4
	m4 $(M4FLAGS) $(ARGS_MANPAGE) $? >$@ || rm $@

# The nfconfig.in file is turned into nfconfig.inc by the configure
# script. The man page is generated on the developers machine for the
# dist.
EXTRA_DIST = nfconfig.in nfconfig.inc $(man_MANS) \
	netcdf_f.h CMakeLists.txt nfconfig.in.cmake

# These files need to be cleaned up my make.
CLEANFILES = nfconfig1.inc netcdf.inc
DISTCLEANFILES = nfconfig.inc
