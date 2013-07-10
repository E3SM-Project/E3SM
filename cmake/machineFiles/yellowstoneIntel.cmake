# CMake initial cache file for yellowstone intel
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpiCC CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PNETCDF} CACHE FILEPATH "")
SET (HDF5 /glade/apps/opt/hdf5/1.8.9/intel/default CACHE FILEPATH "")
SET (SZIP /glade/apps/opt/szlib/2.1/intel/12.1.4 CACHE FILEPATH "")
SET (ZLIB /glade/apps/opt/zlib/1.2.7/intel/12.1.4 CACHE FILEPATH "")
