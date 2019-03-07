# Process HOMME.
SET (HOMME_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/homme)
SET (HOMME_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/homme)
ADD_SUBDIRECTORY(homme)

IF (NOT DEFINED HOMME_TARGET)
  SET (HOMME_TARGET "preqx")
ENDIF()

# Macro to create the lib for dynamics.
# NOTE: You must define PREQX_NP, PREQX_PLEV, PREQX_QSIZE_D, PREQX_USE_PIO, and PREQX_USE_ENERGY
#       BEFORE this point, if you want homme to use custom values (rather than defaults)
MACRO (CreateDynamicsLib)
    
  SET (SCREAM_DYNAMICS_SRC_DIR ${SCREAM_SRC_DIR}/dynamics/homme)

  # It seems redundant to redefine these, when we defined them at the top.
  # However, recall that those variables are local to the folder (and subfolders)
  # where this CMakeLists resides, while this macro could be called from a folder
  # outside this directory.
  SET (HOMME_SOURCE_DIR ${SCREAM_DYNAMICS_SRC_DIR}/homme)
  SET (HOMME_BINARY_DIR ${SCREAM_BINARY_DIR}/src/dynamics/homme/homme)

  SET(SCREAM_DYNAMICS_SOURCES
    ${SCREAM_DYNAMICS_SRC_DIR}/homme_dynamics.cpp
  )

  # CMake variables needed for the config.h file
  SET(NUM_POINTS ${PREQX_NP})
  SET(NUM_PLEV ${PREQX_PLEV})
  SET(NUM_CELLS ${PREQX_NC}) # unused, but needs to be defined
  SET(PIO ${PREQX_USE_PIO})
  SET(PIO_INTERP NOT ${PREQX_USE_PIO})
  SET(ENERGY_DIAGNOSTICS ${PREQX_USE_ENERGY})

  STRING (TOLOWER "${HOMME_TARGET}" HOMME_TARGET)
  IF ("${HOMME_TARGET}" STREQUAL "preqx")
    # Setup the homme target
    PREQX_KOKKOS_SETUP()
  ELSE ()
    # TODO: Add theta when ready in homme
    MESSAGE (FATAL_ERROR "Homme target '${HOMME_TARGET}' not yet supported.")
  ENDIF()

  # This prevents some macros to be defined, such as nested openmp
  SET (USE_KOKKOS_KERNELS TRUE)
  # This informs homme's config files that cmake was used
  # Note: homme already has this in its main CMakeLists, but the scope of ADD_DEFINITIONS
  #       is limited to the directory where it's run and below. This needs to be run
  #       when we create the library, not when we parse homme's directory first
  ADD_DEFINITIONS(-DHOMMEXX_CONFIG_IS_CMAKE)

  # This is needed to compile the library with the correct compile-time options
  SET(THIS_CONFIG_IN ${HOMME_SOURCE_DIR}/src/preqx_kokkos/config.h.cmake.in)
  SET(THIS_CONFIG_HC ${CMAKE_CURRENT_BINARY_DIR}/config.h.c)
  SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)

  # Configure the file
  HommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

  # Append all homme sources to the dynamics sources
  SET(SCREAM_DYNAMICS_SOURCES
      ${SCREAM_DYNAMICS_SOURCES}
      ${PREQX_DEPS_CXX}
      ${PREQX_DEPS_F90}
  )

  # Create the library and set include/link attributes
  ADD_LIBRARY(scream_dynamics ${SCREAM_DYNAMICS_SOURCES})
  target_include_directories(
    scream_dynamics PUBLIC ${SCREAM_SRC_DIR}
    ${CMAKE_CURRENT_BINARY_DIR} # For homme config.h and config.h.c files
    ${EXEC_INCLUDE_DIRS}   # This comes from PREQX_KOKKOS_SETUP
  )

  TARGET_LINK_LIBRARIES(
    scream_dynamics scream_share
    -L${SCREAM_LIBRARY_DIRS} -L${SCREAM_TPL_LIBRARY_DIRS} -L${KOKKOS_LIBRARY_DIR}
    ${SCREAM_TPL_LIBRARIES} ${KOKKOS_TPL_LIBRARIES} ${KOKKOS_LIBRARIES}
    pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} # These come from homme
  )

  SET_TARGET_PROPERTIES(scream_dynamics PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scream_dynamics_modules)
  SET_TARGET_PROPERTIES(scream_dynamics PROPERTIES COMPILE_DEFINITIONS HAVE_CONFIG_H)

ENDMACRO()
