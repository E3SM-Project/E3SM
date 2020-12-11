# CMake initial cache file for Cori KNL nodes

#SET (HOMME_FIND_BLASLAPACK "TRUE" CACHE FILEPATH "")
SET (HOMME_USE_MKL "TRUE" CACHE FILEPATH "") # for Intel

SET (CMAKE_Fortran_COMPILER ftn CACHE FILEPATH "")
SET (CMAKE_C_COMPILER cc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER CC CACHE FILEPATH "")

SET (PIO_FILESYSTEM_HINTS lustre CACHE FILEPATH "")

SET (NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PARALLEL_NETCDF_DIR} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")

SET (MKLROOT $ENV{MKLROOT} CACHE FILEPATH "")

SET (ADD_Fortran_FLAGS "-traceback -craype-verbose" CACHE STRING "")

SET (CMAKE_SYSTEM_NAME Catamount CACHE FILEPATH "")

SET (USE_QUEUING FALSE CACHE BOOL "")

SET (USE_MPIEXEC "srun" CACHE STRING "")

# temporary fix:
SET (USE_MPI_OPTIONS "-c 4 --cpu_bind=cores" CACHE STRING "")

SET (CPRNC_DIR /global/cfs/cdirs/e3sm/tools/cprnc CACHE FILEPATH "")

# by default, cori env loads haswell mod, do
# module unload craype-haswell ; module load craype-mic-knl
#
# NOTE: 2018/2: none of the below seems necessary as long as these modules are loaded:
#  module load cray-netcdf-hdf5parallel
#  module load cray-parallel-netcdf/1.6.1
#
#
# #ZLIB_DIR=/global/common/cori/software/zlib/1.2.8/hsw/intel
# SET (ZLIB_DIR $ENV{ZLIB_DIR} CACHE FILEPATH "")
# SET (ZLIB_LIBRARY ${ZLIB_DIR}/lib/libz.a CACHE FILEPATH "")
#
#EXECUTE_PROCESS(COMMAND which nf-config
#  RESULT_VARIABLE NFCONFIG_PATH_RESULT
#  OUTPUT_VARIABLE NFCONFIG_PATH_OUTPUT
#  ERROR_VARIABLE  NFCONFIG_PATH_ERROR
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#)
#EXECUTE_PROCESS(COMMAND ${NFCONFIG_PATH_OUTPUT} --flibs
#  RESULT_VARIABLE NFCONFIG_RESULT
#  OUTPUT_VARIABLE NFCONFIG_OUTPUT
#  ERROR_VARIABLE  NFCONFIG_ERROR
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#)
#IF (${NFCONFIG_ERROR})
#  MESSAGE(WARNING "${NETCDF_DIR}/bin/nf-config --flibs produced an error. Default linking will be us#ed.")
#ELSE ()
#  SET (ADD_LINKER_FLAGS " ${NFCONFIG_OUTPUT} " CACHE STRING "")
#ENDIF ()
##MESSAGE(STATUS " cori-knl.cmake NFCONFIG_OUTPUT=${NFCONFIG_OUTPUT}")


