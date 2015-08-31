#==============================================================================
#
#  This file sets the environment variables needed to configure and build
#  on the NCAR Wyoming Supercomputing Center systems 
#  (yellowstone/caldera/geyser).
#
#==============================================================================

# Assume all package locations (NetCDF, PnetCDF, HDF5, etc) are already
# set with existing environment variables: NETCDF, PNETCDF, HDF5, etc.

# Define the extra CMake configure options
set (CTEST_CONFIGURE_OPTIONS "-DCMAKE_VERBOSE_MAKEFILE=TRUE")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DPREFER_STATIC=TRUE")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DNetCDF_PATH=$NETCDF_DIR")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DPnetCDF_PATH=$PARALLEL_NETCDF_DIR")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DHDF5_PATH=$HDF5_DIR")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DMPI_C_INCLUDE_PATH=$MPICH_DIR/include")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DMPI_Fortran_INCLUDE_PATH=$MPICH_DIR/include")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DMPI_C_LIBRARIES=$MPICH_DIR/lib/libmpich.a")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DMPI_Fortran_LIBRARIES=$MPICH_DIR/lib/libmpichf90.a")
set (CTEST_CONFIGURE_OPTIONS "${CTEST_CONFIGURE_OPTIONS} -DCMAKE_SYSTEM_NAME=Catamount")
