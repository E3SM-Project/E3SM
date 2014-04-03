# CMake initial cache file for Sandia's climate 48 core linux box
#
#
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

# override cmake's intel defaults:
# default cmake options for Intel: 
#      -assume byterecl -fp-model precise -ftz -g -O3 
# -O3 causes problems on redsky when openMP is enabled (even for 1 thread)
#
#SET (FORCE_Fortran_FLAGS "-traceback -openmp -fp-model fast -ftz -g -O2" CACHE STRING "")
SET (FORCE_Fortran_FLAGS "-traceback -openmp -fp-model precise -ftz -g -O2" CACHE STRING "")
SET (ENABLE_OPENMP TRUE CACHE BOOL "")

# ignored if FORCE_Fortran_FLAGS used above
#SET (ADD_Fortran_FLAGS "-traceback" CACHE STRING "")

SET (NETCDF_DIR "/home/ccsm/netcdf-4.3.0-intel"  CACHE FILEPATH "")
SET (PNETCDF_DIR "/home/ccsm/pnetcdf1.3.1-intel" CACHE FILEPATH "")
#SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

