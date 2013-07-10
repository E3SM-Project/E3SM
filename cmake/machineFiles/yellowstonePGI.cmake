# CMake initial cache file for yellowstone pgi
SET (CMAKE_Fortran_COMPILER mpipf90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpipcc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpipCC CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PNETCDF} CACHE FILEPATH "")
SET (HDF5 /glade/apps/opt/hdf5/1.8.9/pgi/default/ CACHE FILEPATH "")
SET (SZIP /glade/apps/opt/szlib/2.1/pgi/12.5/ CACHE FILEPATH "")
SET (ZLIB /glade/apps/opt/zlib/1.2.7/pgi/12.5/ CACHE FILEPATH "")
