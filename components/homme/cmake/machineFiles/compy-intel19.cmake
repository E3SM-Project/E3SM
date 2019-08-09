# O. Guba 08/2019 
#module purge
#module load cmake/3.11.4
#module load intel/19.0.3
#module load mkl/2019u3
#module load intelmpi/2019u3
#module load netcdf/4.6.3


SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

SET (NETCDF_DIR "$ENV{NETCDF_ROOT}" CACHE STRING "")

#is there a better way? without this linker does not see C netcdf library
SET (ADD_LINKER_FLAGS "-L/$ENV{NETCDF_ROOT}/lib -lnetcdff -lnetcdf" CACHE STRING "")
SET (ADD_Fortran_FLAGS " -traceback " CACHE STRING "")

#SET (USE_MPIEXEC "srun" CACHE STRING "")
#SET (USE_MPI_OPTIONS "--cpu_bind=cores" CACHE STRING "")
SET (HOMME_USE_MKL "TRUE" CACHE FILEPATH "") # for Intel
SET (USE_QUEUING FALSE CACHE BOOL "")
# for standalone HOMME builds:
#SET (CPRNC_DIR /lcrc/group/acme/tools/cprnc CACHE FILEPATH "")

