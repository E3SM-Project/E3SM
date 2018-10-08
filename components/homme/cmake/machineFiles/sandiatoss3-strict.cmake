
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

# Openmpi 1.8 only
SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")

# Openmpi 1.6
#SET (USE_MPI_OPTIONS "-loadbalance" CACHE FILEPATH "")

# This is ignored if we use FORCE_Fortran_FLAGS
# also, these flags will be active since they are going to be last
# like in `ifort: command line warning #10121: overriding '-fp-model fast' with '-fp-model strict'`

SET (ADD_Fortran_FLAGS "-traceback -fp-model strict -qopenmp -O1" CACHE STRING "")
SET (ADD_C_FLAGS "-traceback -fp-model strict -qopenmp -O3" CACHE STRING "")
SET (ADD_CXX_FLAGS "-traceback -fp-model strict -qopenmp -O3" CACHE STRING "")

# redsky upgrade 8/2017, need to load sems-netcdf module:
SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{SEMS_NETCDF_ROOT} CACHE FILEPATH "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET (CPRNC_DIR /projects/ccsm/cprnc/build.toss3 CACHE FILEPATH "")



