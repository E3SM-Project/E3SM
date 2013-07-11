# CMake initial cache file for yellowstone pathscale
SET (CMAKE_Fortran_COMPILER mpief90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpiecc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpiecc CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PNETCDF} CACHE FILEPATH "")
SET (HDF5_DIR /glade/apps/opt/hdf5/1.8.9/pathscale/default CACHE FILEPATH "")
SET (SZIP_DIR /glade/apps/opt/szlib/2.1/pathscale/4.0.12.1/ CACHE FILEPATH "")
SET (ZLIB_DIR /glade/apps/opt/zlib/1.2.7/pathscale/4.0.12.1/ CACHE FILEPATH "")
