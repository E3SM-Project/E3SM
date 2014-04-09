PROJECT(PIO C Fortran)
ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

ENABLE_TESTING()

IF (USER_CMAKE_MODULE_PATH)
  SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${USER_CMAKE_MODULE_PATH})
ELSE()
  SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/pio/cmake")
ENDIF()

find_file( TESTFILE NAMES TryCSizeOf.f90 PATHS ${CMAKE_MODULE_PATH} NO_DEFAULT_PATH)
get_filename_component( TESTFILEPATH ${TESTFILE} PATH)


#IF("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
#  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -wmismatch=mpi_send,mpi_recv,mpi_bcast,mpi_allreduce,mpi_reduce,mpi_isend,mpi_irecv,mpi_irsend,mpi_rsend,mpi_gatherv,mpi_gather,mpi_scatterv,mpi_allgather,mpi_alltoallv,mpi_file_read_all,mpi_file_write_all,mpibcast,mpiscatterv -DCPRNAG")
#ENDIF()

SET(pio_include_dirs_  ${CMAKE_CURRENT_SOURCE_DIR}/pio ${CMAKE_CURRENT_BINARY_DIR}/pio)
SET(PIO_LIB_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/pio ${CMAKE_CURRENT_BINARY_DIR}/pio CACHE STRING "")

#SET(bld_PIO_DEFINITIONS)

TRY_COMPILE(WITH_CSIZEOF ${CMAKE_CURRENT_BINARY_DIR}/tryCompileCSIZEOF 
                        ${TESTFILEPATH}/TryCSizeOf.f90)
#MESSAGE(STATUS "c_sizeof test ${WITH_CSIZEOF}")
IF(${WITH_CSIZEOF} STREQUAL FALSE)
    MESSAGE(STATUS "Fortran compiler does not support c_sizeof function")
    SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -DNO_C_SIZEOF)
endif()


TRY_COMPILE(WITHOUT_MPISERIAL ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPISERIAL
                            ${TESTFILEPATH}/TryMPISERIAL.f90)

IF(${WITHOUT_MPISERIAL} STREQUAL FALSE)
    MESSAGE("Building PIO with mpi serial stubs library")
    SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_MPISERIAL)
ENDIF()


# Netcdf is required

SET (NETCDF_FIND_COMPONENTS F90)
FIND_PACKAGE(NETCDF)
IF (${NETCDF_FOUND})
  MESSAGE("Building PIO with netcdf support ") 
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NETCDF)
  IF (${WITHOUT_MPISERIAL} AND ${NETCDF4_PARALLEL})
    MESSAGE("Building PIO with netcdf4 parallel support") 
    SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NETCDF4)
  ENDIF ()
ENDIF ()
#    SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NETCDF4)



# PNetcdf is optional but used by default
OPTION(WITH_PNETCDF "Whether to build with PNetcdf" TRUE)
IF (${WITH_PNETCDF})
  FIND_PACKAGE(Pnetcdf REQUIRED)
ELSE ()
  MESSAGE(WARNING "Warning: Not building with PNetcdf - cannot run all regression tests.")
ENDIF ()


OPTION(PIO_BIG_ENDIAN "Specify that the machine is big endian" test)
IF ("${PIO_BIG_ENDIAN}" STREQUAL "test") 
  INCLUDE(TestBigEndian)
  TestBigEndian(PIO_BIG_ENDIAN_TEST)
  IF(PIO_BIG_ENDIAN_TEST) 
    SET(PIO_BIG_ENDIAN ON  CACHE BOOL "")
  ELSE()
    SET(PIO_BIG_ENDIAN OFF CACHE BOOL "")
  ENDIF()
ELSE()
  SET(PIO_BIG_ENDIAN ${PIO_BIG_ENDIAN} CACHE BOOL "")
ENDIF()                                            


IF (PIO_FILESYSTEM_HINTS STREQUAL "lustre")
  MESSAGE(STATUS "PIO using lustre filesystem hints")
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -DPIO_LUSTRE_HINTS)
ELSEIF(PIO_FILESYSTEM_HINTS STREQUAL "gpfs")
  MESSAGE(STATUS "PIO using gpfs filesystem hints")
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -DPIO_GPFS_HINTS)
ELSEIF(NOT "${PIO_FILESYSTEM_HINTS}" STREQUAL "")
  MESSAGE(WARNING "${PIO_FILESYSTEM_HINTS} not valid option for PIO_FILESYSTEM_HINTS; use gpfs or lustre.")
ENDIF()

IF(NETCDF_FOUND) 
  SET(pio_include_dirs_ ${pio_include_dirs_} ${NETCDF_INCLUDE_DIR})
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NETCDF)
ELSE()
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NONETCDF)
ENDIF()
IF(PNETCDF_FOUND)
  SET(pio_include_dirs_ ${pio_include_dirs_} ${PNETCDF_INCLUDE_DIR})
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_PNETCDF)
ELSE()
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NOPNETCDF)
ENDIF()

FIND_PACKAGE(MPI)

OPTION(PIO_USE_MPIIO "Enable support for MPI-IO (default auto detect)" ON)
IF (PIO_USE_MPIIO) 
  TRY_COMPILE(MPIIO_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPIIO 
                            ${TESTFILEPATH}/TryMPIIO.f90)
  IF (${MPIIO_SUCCESS})
    MESSAGE(STATUS "MPIIO detected and enabled.")
  ELSE()
    MESSAGE(STATUS "MPIIO not detected and therefore disabled.")
    SET(PIO_USE_MPIIO FALSE)
  ENDIF()
ENDIF()
IF (${PIO_USE_MPIIO}) 
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -DUSEMPIIO)
ENDIF()

SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NOUSEMCT)


OPTION(PIO_USE_MPIMOD "Use Fortran MPI module (default auto detect)" ON)
IF (PIO_USE_MPIMOD)
  TRY_COMPILE(MPIMOD_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPIMod
                             ${TESTFILEPATH}/TryMPIMod.f90)
  IF (${MPIMOD_SUCCESS})
    MESSAGE(STATUS "MPI Fortran module detected and enabled.")
  ELSE()
    MESSAGE(STATUS "MPI Fortran module not detected and therefore disabled.")
    SET(PIO_USE_MPIMOD FALSE)
  ENDIF()
ENDIF()

IF (NOT ${PIO_USE_MPIMOD})
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -DNO_MPIMOD)
ENDIF()


OPTION(PIO_BUILD_TIMING  "" OFF)

if(${PIO_BUILD_TIMING})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/timing)
    SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -DTIMING -I ${CMAKE_CURRENT_BINARY_DIR}/timing)
    ADD_SUBDIRECTORY(timing timing)
  endif()
endif()

SET(PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} CACHE STRING "")
ADD_DEFINITIONS(${PIO_DEFINITIONS})

SET(PIO_INCLUDE_DIRS   ${pio_include_dirs_} CACHE STRING "")
INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})

ADD_SUBDIRECTORY(pio)
ADD_SUBDIRECTORY(test)
