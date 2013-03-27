PROJECT(HOMME C CXX Fortran)

ENABLE_LANGUAGE(Fortran) 
INCLUDE(FortranCInterface)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# CTest For Nightly Builds
OPTION(HOMME_USE_CTEST FALSE)

IF (${HOMME_USE_CTEST})
  INCLUDE(CTest)
ENDIF ()

# For testing
ENABLE_TESTING()

# Add the cmake/CMakeLists.txt file
SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


# need to get the cmake_fortran_c_interface.h header
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
FortranCInterface_HEADER(cmake_fortran_c_interface.h
                         MACRO_NAMESPACE  "FCI_")

# Set up some directories to facilitate finding netcdf and pnetcdf
INCLUDE(SetHostPaths)
INCLUDE(HommeMacros)

# Default don't use quad precision
OPTION(HOMME_QUAD_PREC "" TRUE)

IF (${HOMME_QUAD_PREC}) 
  testQuadPrec(HOMME_QUAD_SUPPORTED)
  IF (NOT ${HOMME_QUAD_SUPPORTED})
    SET(HOMME_QUAD_PREC FALSE)
    MESSAGE(STATUS "Quadruple-precision requested but unavailable on this
                    system with this compiler")
  ENDIF()
ENDIF ()

#MESSAGE(STATUS "CMAKE_IMP_LIBS=${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")

# CPRNC for diffing Netcdf files
OPTION(TEST_USING_CPRNC "Whether to use cprnc for testing" FALSE)

IF (${TEST_USING_CPRNC})
  find_package(Cprnc REQUIRED)
  MESSAGE(STATUS "Using CPRNC to diff the Netcdf files")
ELSE ()
  MESSAGE(STATUS 
    "Not using CPRNC to diff the Netcdf files, only diffing stdout")
ENDIF ()

# Use python for 
find_package(PythonInterp REQUIRED)

# need to tell the source to use this capability
ADD_DEFINITIONS(-DINCLUDE_CMAKE_FCI)

# Some options for conditional builds
OPTION(BUILD_HOMME_SWEQX "" ON)
OPTION(BUILD_HOMME_PREQX "" ON)
OPTION(BUILD_HOMME_SWDGX "" ON)
OPTION(BUILD_HOMME_SWIM  "" OFF)

IF (${SWEQX_ONLY})
  SET(BUILD_HOMME_SWEQX "" ON)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
ENDIF ()

IF (${PREQX_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" ON)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" OFF)
ENDIF ()

IF (${SWDGX_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_SWDGX "" ON)
  SET(BUILD_HOMME_SWIM  "" OFF)
ENDIF ()

IF (${SWIM_ONLY})
  SET(BUILD_HOMME_SWEQX "" OFF)
  SET(BUILD_HOMME_PREQX "" OFF)
  SET(BUILD_HOMME_SWDGX "" OFF)
  SET(BUILD_HOMME_SWIM  "" ON)
ENDIF ()

# to allow ifdefs to be included
if(REFSOLN)
  MESSAGE(STATUS "Including option to read or write reference soln data...")
  ADD_DEFINITIONS(-D_REFSOLN)
endif()

#ADD_SUBDIRECTORY(cmake)

MESSAGE(STATUS "Homme hostname = ${Homme_Hostname}")
MESSAGE(STATUS "Homme OS = ${Homme_OS}")
MESSAGE(STATUS "C Compiler= ${CMAKE_C_COMPILER_ID}")
MESSAGE(STATUS "CXX Compiler= ${CMAKE_CXX_COMPILER_ID}")
MESSAGE(STATUS "Fortran Compiler= ${CMAKE_Fortran_COMPILER_ID}")
MESSAGE(STATUS "CMAKE_SYSTEM= ${CMAKE_SYSTEM}")
MESSAGE(STATUS "CMAKE_SYSTEM_NAME= ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "HOMME_FIND_BLASLAPACK=${HOMME_FIND_BLASLAPACK}")

OPTION(HOMME_ENABLE_OPENMP FALSE)

IF (HOMME_ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    # The fortran openmp flag should be the same as the C Flag
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  ELSE ()
    MESSAGE(FATAL_ERROR "Unable to find OpenMP")
  ENDIF()
ENDIF ()

# For PIO Set these options to allow a longer line length (until PIO is fixed)
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
ELSEIF (CMAKE_COMPILER_ID STREQUAL PGI)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mextend")
ELSEIF (CMAKE_COMPILER_ID STREQUAL PathScale)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend-source")
ENDIF ()

# Intel specific flags
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise -ftz")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O3")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2")
ENDIF ()

IF(HOMME_FIND_BLASLAPACK)
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ELSE()
  SET(HOMME_FIND_BLASLAPACK FALSE)
  MESSAGE(STATUS "HOMME_FIND_BLASLAPACK=${HOMME_FIND_BLASLAPACK}")
  # Build the blas and lapack provided with homme
  ADD_SUBDIRECTORY(libs/blas)
  ADD_SUBDIRECTORY(libs/lapack)
ENDIF()

FIND_PACKAGE(Netcdf REQUIRED)
IF (${WITH_PNETCDF})
  FIND_PACKAGE(Pnetcdf REQUIRED)
ELSE ()
  FIND_PACKAGE(Pnetcdf)
ENDIF ()

INCLUDE(ExternalProject)

get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(CACHE_VAR ${CACHE_VARS})
  #get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
#  if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
      set(CACHE_VAR_TYPE)
    else()
      set(CACHE_VAR_TYPE :${CACHE_VAR_TYPE})
    endif()
    set(HOMME_CMAKE_ARGS ${HOMME_CMAKE_ARGS} "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}")
  #endif()
endforeach()
message("HOMME_CMAKE_ARGS: ${HOMME_CMAKE_ARGS}")

# CPRNC
SET (CPRNC_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/utils/cprnc)
SET (CPRNC_BINARY ${CMAKE_CURRENT_BINARY_DIR}/utils/cprnc/bin/cprnc)

# Compiler options
#SET(COMPILE_ARGS ${COMPILE_ARGS} 
#    "-DCMAKE_Fortran_COMPILER:FILEPATH=\"${CMAKE_Fortran_COMPILER}\""
#    "-DCMAKE_C_COMPILER:FILEPATH=\"${CMAKE_C_COMPILER}\""
#    "-DCMAKE_CXX_COMPILER:FILEPATH=\"${CMAKE_CXX_COMPILER}\"")
SET(COMP_FORT "-DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}")
SET(COMP_CC "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}")

SET(COMP_CXX "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}")


string (REPLACE ";" " " COMPILE_ARGS "${COMPILE_ARGS}")
#IF (DEFINED CPRNC_NO_RECURSE AND ${CPRNC_NO_RECURSE})
# The following builds cprnc 
ExternalProject_Add(
  cprnc
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils/cprnc
  CMAKE_ARGS
    -DNETCDF_DIR=${CPRNC_NETCDF_DIR}
    -DHDF5_DIR=${CPRNC_HDF5_DIR}
    -DSZIP_DIR=${CPRNC_SZIP_DIR}
    -DZLIB_DIR=${CPRNC_ZLIB_DIR}
    -DCURL_DIR=${CPRNC_CURL_DIR}
    -DCMAKE_INSTALL_PREFIX=${CPRNC_INSTALL_DIR}
  INSTALL_DIR ${CPRNC_INSTALL_DIR}
  LOG_CONFIGURE ON)
#ENDIF ()

#SET (HOMME_TAG_DIR ${CMAKE_CURRENT_BINARY_DIR}/homme_tag)


##MESSAGE(STATUS "HOMME_CMAKE_ARGS=${HOMME_CMAKE_ARGS}")
##string (REPLACE ";" " " HOMME_CMAKE_ARGS "${HOMME_CMAKE_ARGS}")

#MESSAGE(STATUS "HOMME_CMAKE_ARGS=${HOMME_CMAKE_ARGS}")
#MESSAGE(STATUS "COMPILE_ARGS = ${COMPILE_ARGS}")

#ExternalProject_Add(
#    homme_tag
#    SVN_REPOSITORY https://svn-homme-model.cgd.ucar.edu/trunk
#    CMAKE_ARGS 
#      ${HOMME_CMAKE_ARGS}
#      -DCMAKE_INSTALL_PREFIX=${HOMME_TAG_DIR}
#    INSTALL_DIR ${HOMME_TAG_DIR}
#    LOG_CONFIGURE ON
#    CMAKE_CACHE_ARGS 
#      ${COMP_FORT}
#      ${COMP_CC}
#      ${COMP_CXX}
#    LOG_DOWNLOAD ON
#    LOG_CONFIGURE ON
#    LOG_BUILD ON)
#ExternalProject_Add(
#    homme_tag
#    SVN_REPOSITORY https://svn-homme-model.cgd.ucar.edu/trunk
#    CMAKE_ARGS 
#      -DCMAKE_INSTALL_PREFIX=${HOMME_TAG_DIR}
#    INSTALL_DIR ${HOMME_TAG_DIR}
#    LOG_CONFIGURE ON
#    CMAKE_CACHE_ARGS 
#      ${HOMME_CMAKE_ARGS}
#    LOG_DOWNLOAD ON
#    LOG_CONFIGURE ON
#    LOG_BUILD ON)


ADD_SUBDIRECTORY(utils/pio)
ADD_SUBDIRECTORY(utils/timing)

if(${BUILD_HOMME_SWEQX})
  ADD_SUBDIRECTORY(src/sweqx)
endif()

if(${BUILD_HOMME_PREQX})
  ADD_SUBDIRECTORY(src/preqx)
endif()

if(${BUILD_HOMME_SWDGX})
  ADD_SUBDIRECTORY(src/swdgx)
endif()

if(${BUILD_HOMME_SWIM})
  MESSAGE(STATUS "Building swim...")
  FIND_PACKAGE(Trilinos REQUIRED)
  MESSAGE(STATUS "Using Trilinos: ${Trilinos_LIBRARY_DIRS}")
  ADD_SUBDIRECTORY(src/swim)
endif()

ADD_SUBDIRECTORY(test_execs)
