# CMake initial cache file for titan
SET (CMAKE_Fortran_COMPILER mpixlf90_r CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpixlc_r CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpixlcxx_r CACHE FILEPATH "")
SET (NETCDF_DIR /u/an/es/bjamroz/software/netcdf-4.2 CACHE FILEPATH "")
SET (PNETCDF_DIR /u/an/es/bjamroz/software/parallel-netcdf-1.4.1 CACHE FILEPATH "")

SET (USE_MPIEXEC srun CACHE FILEPATH "")

SET (IS_BIG_ENDIAN TRUE CACHE BOOL "")
