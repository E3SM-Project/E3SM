# O.Guba 08/2019
# module load cmake intel/18.0.0 intelmpi/2018 mkl netcdf
# however, NETCDF_ROOT won't be set

SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

#intel18
SET (NETCDF_DIR "/share/apps/netcdf/4.6.3/intel/18.0.0/" CACHE STRING "")

#duplicated path
#is there a better way? without this linker does not see C netcdf library
SET (ADD_LINKER_FLAGS "-L//share/apps/netcdf/4.6.3/intel/18.0.0//lib -lnetcdff -lnetcdf" CACHE STRING "")
#SET (ADD_Fortran_FLAGS " -traceback " CACHE STRING "")

SET (USE_MPIEXEC "srun" CACHE STRING "")
SET (USE_MPI_OPTIONS " --mpi=pmi2 --kill-on-bad-exit --cpu_bind=cores" CACHE STRING "")
SET (HOMME_USE_MKL "TRUE" CACHE FILEPATH "") # for Intel
SET (USE_QUEUING FALSE CACHE BOOL "")
# for standalone HOMME builds:
#SET (CPRNC_DIR /lcrc/group/acme/tools/cprnc CACHE FILEPATH "")

