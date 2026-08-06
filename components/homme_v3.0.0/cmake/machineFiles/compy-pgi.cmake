# O. Guba 08/2019 
#module purge

SET(CMAKE_C_COMPILER "mpipgcc" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpipgcxx" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpipgf90" CACHE STRING "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

SET (NETCDF_DIR "$ENV{NETCDF_ROOT}" CACHE STRING "")

#is there a better way? without this linker does not see C netcdf library
SET (ADD_LINKER_FLAGS "-L/$ENV{MKLROOT}/lib/intel64 -L/$ENV{NETCDF_ROOT}/lib -lnetcdff -lnetcdf -lmkl_rt" CACHE STRING "")
SET (ADD_Fortran_FLAGS "-traceback" CACHE STRING "")

SET (USE_MPIEXEC "srun" CACHE STRING "")
SET (USE_MPI_OPTIONS "--mpi=pmi2 --kill-on-bad-exit --cpu_bind=cores" CACHE STRING "")
SET (USE_QUEUING FALSE CACHE BOOL "")
# for standalone HOMME builds:
#SET (CPRNC_DIR /lcrc/group/acme/tools/cprnc CACHE FILEPATH "")
