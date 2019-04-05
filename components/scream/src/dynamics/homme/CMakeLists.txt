# Process HOMME.
# NOTE: You must define PREQX_NP, PREQX_PLEV, PREQX_QSIZE_D, PREQX_USE_PIO, and PREQX_USE_ENERGY
#       BEFORE this point, if you want homme to use custom values (rather than defaults)

SET (HOMME_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/homme)
SET (HOMME_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/homme)

# Disable all the targets by default
SET(BUILD_HOMME_SWEQX        OFF CACHE BOOL "")
SET(BUILD_HOMME_PREQX        OFF CACHE BOOL "")
SET(BUILD_HOMME_THETA        OFF CACHE BOOL "")
SET(BUILD_HOMME_PREQX_ACC    OFF CACHE BOOL "")
SET(BUILD_HOMME_PREQX_KOKKOS OFF CACHE BOOL "")
SET(BUILD_HOMME_PESE         OFF CACHE BOOL "")
SET(BUILD_HOMME_SWIM         OFF CACHE BOOL "")
SET(BUILD_HOMME_PRIM         OFF CACHE BOOL "")

# Default homme target is preqx, for now
IF (NOT DEFINED HOMME_TARGET)
  SET (HOMME_TARGET "preqx_kokkos")
ENDIF()
STRING (TOLOWER "${HOMME_TARGET}" HOMME_TARGET)

# Enable the requested target
IF ("${HOMME_TARGET}" STREQUAL "preqx_kokkos")
  SET (BUILD_HOMME_PREQX_KOKKOS TRUE)
ELSE ()
  MESSAGE (FATAL_ERROR "Homme target '${HOMME_TARGET}' is not supported within scream.\n")
ENDIF()

ADD_SUBDIRECTORY(homme)

MACRO (CreateDynamicsLib libName)
  ADD_DEFINITIONS(-DHOMMEXX_CONFIG_IS_CMAKE)
  SET (SCREAM_DYNAMICS_SRC_DIR ${SCREAM_SRC_DIR}/dynamics/homme)

  # Re-set these, so that they are visible from wherever the macro is called
  SET (HOMME_SOURCE_DIR ${SCREAM_SOURCE_DIR}/src/dynamics/homme/homme)
  SET (HOMME_BINARY_DIR ${SCREAM_BINARY_DIR}/src/dynamics/homme/homme)

  IF ("${HOMME_TARGET}" STREQUAL "preqx_kokkos")
    # Setup the homme target
    PREQX_KOKKOS_SETUP()
  ELSE ()
    MESSAGE (FATAL_ERROR "Homme target '${HOMME_TARGET}' is not supported within scream.\n \
                          Did you forget to set 'HOMME_TARGET' before calling 'CreateDynamicsLib'?\n")
  ENDIF()

  # Get homme include dirs. This var is set by the ${HOMME_TARGET}_setup macro
  SET (HOMME_INCLUDE_DIRS ${EXEC_INCLUDE_DIRS})

  SET (SCREAM_DYNAMICS_SOURCES
       ${SCREAM_DYNAMICS_SRC_DIR}/atmosphere_dynamics.cpp
       ${SCREAM_DYNAMICS_SRC_DIR}/homme_dynamics.cpp
       # ${SCREAM_DYNAMICS_SRC_DIR}/scream_homme_interface.F90
       ${PREQX_DEPS_CXX}
       ${PREQX_DEPS_F90}
  )

  # Append all libraries need byhomme to the list of libraries needed by dynamics
  SET (SCREAM_DYNAMICS_LINK_LIBRARIES
       scream_share
       pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
  )

  # Set the variables that need to be expanded by cmake during configure_file
  SET (NUM_POINTS ${PREQX_NP})
  SET (NUM_CELLS ${PREQX_NC})
  SET (NUM_PLEV ${PREQX_PLEV})
  SET (QSIZE_D ${PREQX_QSIZE_D})

  SET (PIO ${PREQX_USE_PIO})
  SET (PIO_INTERP NOT ${PREQX_USE_PIO})
  SET (ENERGY_DIAGNOSTICS ${PREQX_USE_ENERGY})

  # This is needed to compile the test executables with the correct options for homme
  SET(THIS_CONFIG_IN ${HOMME_SOURCE_DIR}/src/${HOMME_TARGET}/config.h.cmake.in)
  SET(THIS_CONFIG_HC ${CMAKE_CURRENT_BINARY_DIR}/config.h.c)
  SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)

  # First configure the file (which formats the file as C)
  HommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

  ADD_DEFINITIONS(-DHAVE_CONFIG_H)

  ADD_LIBRARY(${libName} ${SCREAM_DYNAMICS_SOURCES})
  SET_TARGET_PROPERTIES (${libName} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${libName}_modules)
  TARGET_INCLUDE_DIRECTORIES(
      ${libName} PUBLIC
      ${SCREAM_INCLUDE_DIRS}
      ${SCREAM_TPL_INCLUDE_DIRS}
      ${EXEC_INCLUDE_DIRS}        # This comes from PREQX_KOKKOS_SETUP
      ${CMAKE_CURRENT_BINARY_DIR} # For config.h
  )

ENDMACRO(CreateDynamicsLib)

ADD_SUBDIRECTORY (tests)
