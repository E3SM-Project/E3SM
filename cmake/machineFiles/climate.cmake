# CMake initial cache file for Sandia's climate 48 core linux box
# 
# disable openMP support
# code compiles and runs with OPENMP support, 1 thread, but
# hangs if using >1 threads.  
#
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

# override cmake's intel defaults:
# default cmake options for Intel: 
#      -assume byterecl -fp-model precise -ftz -g -O3 
# -O3 causes problems on redsky when openMP is enabled (even for 1 thread)
#
#SET (FORCE_Fortran_FLAGS "-openmp -fp-model fast -ftz -g -O2" CACHE STRING "")
#SET (ENABLE_OPENMP TRUE CACHE BOOL "")

SET (ADD_Fortran_FLAGS "-traceback" CACHE STRING "")

SET (NETCDF_DIR $ENV{NETCDF_PATH} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PNETCDF_PATH} CACHE FILEPATH "")
SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

