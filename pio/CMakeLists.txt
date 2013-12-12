IF( NOT GENF90_PATH) 
  SET (GENF90_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
ENDIF()

PROJECT(PIO C Fortran)
ENABLE_LANGUAGE(Fortran)
#INCLUDE(FortranCInterface)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
IF (USER_CMAKE_MODULE_PATH)
  SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${USER_CMAKE_MODULE_PATH})
ELSE()
  SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
ENDIF()


SET(pio_include_dirs_  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
SET(PIO_LIB_DIR  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} CACHE STRING "")

#SET(bld_PIO_DEFINITIONS)

TRY_COMPILE(WITH_CSIZEOF ${CMAKE_CURRENT_BINARY_DIR}/tryCompileCSIZEOF
                            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryCSizeOf.f90)
#MESSAGE(STATUS "c_sizeof test ${WITH_CSIZEOF}")
IF(${WITH_CSIZEOF} STREQUAL FALSE)
    MESSAGE(STATUS "Fortran compiler does not support c_sizeof function")
    SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -DNO_C_SIZEOF)
endif()


TRY_COMPILE(WITHOUT_MPISERIAL ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPISERIAL
                            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryMPISERIAL.f90)

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

OPTION(PIO_USE_MPIIO "Enable support for MPI-IO (default auto detect)" ON)
IF (PIO_USE_MPIIO) 
  TRY_COMPILE(MPIIO_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPIIO 
                            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryMPIIO.f90)
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

OPTION(PIO_USE_BOX  "" ON)
if(PIO_USE_BOX)
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_USEBOX)
else()
  SET(bld_PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} -D_NOUSEBOX)
endif()

OPTION(PIO_USE_MPIMOD "Use Fortran MPI module (default auto detect)" ON)
IF (PIO_USE_MPIMOD)
  TRY_COMPILE(MPIMOD_SUCCESS ${CMAKE_CURRENT_BINARY_DIR}/tryCompileMPIMod
                             ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryMPIMod.f90)
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

SET(PIO_DEFINITIONS ${bld_PIO_DEFINITIONS} CACHE STRING "")
ADD_DEFINITIONS(${PIO_DEFINITIONS})

SET(PIO_INCLUDE_DIRS   ${pio_include_dirs_} CACHE STRING "")
INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})

SET(SRCS_C topology.c)

SET(SRCS_F90 pio.F90 pio_kinds.F90 nf_mod.F90  ionf_mod.F90 pio_types.F90
             piolib_mod.F90 pio_mpi_utils.F90 pio_nf_utils.F90 pio_utils.F90
             pio_support.F90                                calcdisplace_mod.F90
	     calcdecomp.F90 pio_msg_mod.F90 pio_msg_callbacks.F90)

SET(TEMPSRCF90    pionfatt_mod.F90
         	                pionfread_mod.F90
                                pionfwrite_mod.F90
                                pionfput_mod.F90
                                pionfget_mod.F90
	                        alloc_mod.F90
                                box_rearrange.F90
                                rearrange.F90
         	                iompi_mod.F90
           	                piodarray.F90
	                        pio_spmd_utils.F90
	                        pio_msg_getput_callbacks.F90
        )

FOREACH(tempfile IN LISTS TEMPSRCF90)
ADD_CUSTOM_COMMAND(
	OUTPUT ${tempfile}
	COMMAND ${GENF90_PATH}/genf90.pl ${CMAKE_CURRENT_SOURCE_DIR}/${tempfile}.in > ${tempfile}
	DEPENDS  ${tempfile}.in
)
ENDFOREACH()
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  SET(CMAKE_Fortran_FLAGS	"${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif()

ADD_LIBRARY(pio ${SRCS_F90} ${SRCS_C} ${TEMPSRCF90})
if(PNETCDF_LIBRARY)
  TARGET_LINK_LIBRARIES(pio ${PNETCDF_LIBRARIES})
ENDIF()
if(NETCDF_LIBRARIES)
  TARGET_LINK_LIBRARIES(pio ${NETCDF_LIBRARIES})
ENDIF()
TARGET_LINK_LIBRARIES(pio ${ADDITIONAL_LIBS})

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../unittests)
  ADD_SUBDIRECTORY(../unittests unittests)
endif()
#if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../testpio)
#  ADD_SUBDIRECTORY(../testpio testpio)
#endif()
