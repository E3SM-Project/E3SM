# CMake initial cache file for yellowstone pathscale
SET (CMAKE_Fortran_COMPILER mpief90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpiecc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpiecc CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PNETCDF} CACHE FILEPATH "")
SET (HDF5 /glade/apps/opt/hdf5/1.8.9/pathscale/default CACHE FILEPATH "")
SET (SZIP /glade/apps/opt/szlib/2.1/pathscale/4.0.12.1/ CACHE FILEPATH "")
SET (ZLIB /glade/apps/opt/zlib/1.2.7/pathscale/4.0.12.1/ CACHE FILEPATH "")
