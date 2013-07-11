# CMake initial cache file for yellowstone intel
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpiCC CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF})
SET (PNETCDF_DIR $ENV{PNETCDF})
SET (HDF5_DIR /glade/apps/opt/hdf5/1.8.9/intel/default)
SET (SZIP_DIR /glade/apps/opt/szlib/2.1/intel/12.1.4)
SET (ZLIB_DIR /glade/apps/opt/zlib/1.2.7/intel/12.1.4)
