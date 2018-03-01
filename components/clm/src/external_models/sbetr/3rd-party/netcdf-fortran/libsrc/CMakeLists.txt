CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Tell make to install these as a header file.
SET(include_HEADERS nfconfig.inc)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# The (temporary) library we are building.
SET(ncfortran_LIBRARY_SOURCES
  fort-attio.c fort-control.c fort-dim.c
  fort-genatt.c fort-geninq.c fort-genvar.c fort-lib.c fort-misc.c
  fort-v2compat.c fort-vario.c fort-var1io.c fort-varaio.c fort-varmio.c
  fort-varsio.c cfortran.h fort-lib.h ncfortran.h)
# Is the user building netCDF-4?
#IF(USE_NETCDF4)
#  SET(ncfortran_LIBRARY_SOURCES ${ncfortran_LIBRARY_SOURCES} fort-nc4.c)
#ENDIF(USE_NETCDF4)

# The F77 API man page.
SET(man_MANS netcdf_f77.3)

# This rule generates the F77 manpage.
SET(ARGS_MANPAGE "-DAPI=FORTRAN")
IF (USE_NETCDF4)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} "-DNETCDF4=TRUE")
ENDIF(USE_NETCDF4)
IF(BUILD_DAP)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} "-DDAP=TRUE")
ENDIF(BUILD_DAP)
IF(BUILD_PARALLEL)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} "-DPARALLEL_IO=TRUE")
ENDIF(BUILD_PARALLEL)

# TODO: Figure out the cmake syntax for this
# netcdf_f77.3: $(top_srcdir)/docs/netcdf.m4
#	m4 $(M4FLAGS) $(ARGS_MANPAGE) $? >$@ || rm $@

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${CMAKE_CPP_FLAGS})

ADD_LIBRARY(ncfortran ${ncfortran_LIBRARY_SOURCES})
TARGET_LINK_LIBRARIES(ncfortran ${NETCDF_C_LIBRARY} ${EXTRA_DEPS})
# The nfconfig.in file is turned into nfconfig.inc by the configure
# script. The man page is generated on the developers machine for the
# dist.
SET(EXTRA_DIST nfconfig.in nfconfig.inc ${man_MANS} netcdf_f.h)# These files need to be cleaned up my make.

# These files need to be cleaned up my make.
SET(CLEANFILES nfconfig1.inc netcdf.inc)
SET(DISTCLEANFILES nfconfig.inc )
