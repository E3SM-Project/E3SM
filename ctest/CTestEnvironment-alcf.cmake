#==============================================================================
#
#  This file sets the environment variables needed to configure and build
#  on the Argonne Leadership Computing Facility systems 
#  (mira/cetus/vesta/cooley).
#
#==============================================================================

# Assume all package locations (NetCDF, PnetCDF, HDF5, etc) are already
# set with existing environment variables: NETCDF, PNETCDF, HDF5, etc.

# Define the extra CMake configure options
set (CTEST_CONFIGURE_OPTIONS "-DCMAKE_VERBOSE_MAKEFILE=TRUE")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DPREFER_STATIC=TRUE")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DLIBZ_PATH=/soft/libraries/alcf/current/xl/ZLIB/")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DHDF5_PATH=/soft/libraries/hdf5/current/cnk-xl/current")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DNetCDF_PATH=/soft/libraries/netcdf/current/cnk-xl/current")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DPnetCDF_PATH=/soft/libraries/pnetcdf/current/cnk-xl/current")