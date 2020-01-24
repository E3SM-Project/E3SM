# Process HOMME.
# NOTE: You must define PREQX_NP, PREQX_PLEV, PREQX_QSIZE_D, PREQX_USE_PIO, and PREQX_USE_ENERGY
#       BEFORE this point, if you want homme to use custom values (rather than defaults)

SET (HOMME_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/homme CACHE INTERNAL "Homme source directory")
SET (HOMME_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/homme CACHE INTERNAL "Homme binary directory")

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
SET (HOMME_TARGETS "preqx_kokkos;theta_kokkos")

# Enable the requested target
FOREACH (target IN ITEMS ${HOMME_TARGETS})
  IF ("${target}" STREQUAL "preqx_kokkos")
    SET (BUILD_HOMME_PREQX_KOKKOS TRUE)
  ELSEIF ("${target}" STREQUAL "theta_kokkos")
    SET (BUILD_HOMME_THETA_KOKKOS TRUE)
  ELSE ()
    MESSAGE (FATAL_ERROR "Homme target '${target}' is not supported within scream.\n")
  ENDIF()
ENDFOREACH()

SET (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)
ADD_SUBDIRECTORY(homme)

# Take the variable that homme set, and mark them as (internal) cached, so they are globally visible.
# This is needed since the macro CreateDynamicsLib will be invoked from outside this directory,
# for instance, from the scream/tests folder.
SET (HOMME_C_FLAGS       ${HOMME_C_FLAGS}       CACHE INTERNAL "List of C compiler flags for Homme")
SET (HOMME_CXX_FLAGS     ${HOMME_CXX_FLAGS}     CACHE INTERNAL "List of CXX compiler flags for Homme")
SET (HOMME_Fortran_FLAGS ${HOMME_Fortran_FLAGS} CACHE INTERNAL "List of Fortran compiler flags for Homme")

SET (DynamicsLibsCreated   "" CACHE INTERNAL "List of created dynamics libraries names")

MACRO (CreateDynamicsLib HOMME_TARGET NP PLEV QSIZE USE_PIO USE_ENERGY)
  SET (PREQX_NP ${NP})
  SET (PREQX_PLEV ${PLEV})
  SET (PREQX_QSIZE_D ${QSIZE})
  SET (PREQX_USE_PIO ${USE_PIO})
  SET (PREQX_USE_ENERGY ${USE_ENERGY})

  SET (dynLibName "${HOMME_TARGET}_np${NP}_nlev${PLEV}_qsize${QSIZE}_pio${USE_PIO}_energy${USE_ENERGY}")
  IF ("${dynLibName}" IN_LIST DynamicsLibsCreated)
    SET (dynLibDir ${CMAKE_BINARY_DIR}/src/dynamics/homme/${dynLibName})
  ELSE ()
    # Add this library name to the list of created libraries
    SET (DynamicsLibsCreated
         ${DynamicsLibsCreated}
         ${dynLibName}
         CACHE INTERNAL "List created dynamics libraries names"
    )

    # Set the dynamics library directory
    SET (dynLibDir ${CMAKE_BINARY_DIR}/src/dynamics/homme/${dynLibName})

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
         ${SCREAM_DYNAMICS_SRC_DIR}/dynamics_driven_grids_manager.cpp
         ${SCREAM_DYNAMICS_SRC_DIR}/scream_homme_interface.F90
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
    SET(THIS_CONFIG_HC ${dynLibDir}/config.h.c)
    SET(THIS_CONFIG_H ${dynLibDir}/config.h)

    # First configure the file (which formats the file as C)
    HommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

    ADD_DEFINITIONS(-DHAVE_CONFIG_H)

    # Convert strings to lists
    # Note: TARGET_COMPILE_OPTIONS perform 'de-duplication' of repeated flags.
    #       this causes the repeated "-fp-model X -fp-model Y" to be turned into
    #       "-fp-model X Y", which the compiler does not like.
    #       Thankfully, there is a remedy, albeit a bit ugly: prepend the flag to be kept
    #       (even if repeated) with 'SHELL: '.
    STRING(REPLACE "-fp-model" "SHELL: -fp-model" HOMME_C_FLAGS_LIST_TMP       "${HOMME_C_FLAGS}")
    STRING(REPLACE "-fp-model" "SHELL: -fp-model" HOMME_CXX_FLAGS_LIST_TMP     "${HOMME_CXX_FLAGS}")
    STRING(REPLACE "-fp-model" "SHELL: -fp-model" HOMME_Fortran_FLAGS_LIST_TMP "${HOMME_Fortran_FLAGS}")
    STRING(REPLACE " " ";" HOMME_C_FLAGS_LIST       "${HOMME_C_FLAGS_TMP}")
    STRING(REPLACE " " ";" HOMME_CXX_FLAGS_LIST     "${HOMME_CXX_FLAGS_TMP}")
    STRING(REPLACE " " ";" HOMME_Fortran_FLAGS_LIST "${HOMME_Fortran_FLAGS_TMP}")

    ADD_LIBRARY(${dynLibName} ${SCREAM_DYNAMICS_SOURCES})
    target_link_libraries(${dynLibName} ${SCREAM_DYNAMICS_LINK_LIBRARIES})
    SET_TARGET_PROPERTIES (${dynLibName} PROPERTIES Fortran_MODULE_DIRECTORY ${dynLibDir}/${dynLibName}_modules)

    # Add flags, using a generator to achieve language specific behavior
    TARGET_COMPILE_OPTIONS(${dynLibName} BEFORE PUBLIC $<$<COMPILE_LANGUAGE:C>:${HOMME_C_FLAGS_LIST}>)
    TARGET_COMPILE_OPTIONS(${dynLibName} BEFORE PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${HOMME_CXX_FLAGS_LIST}>)
    TARGET_COMPILE_OPTIONS(${dynLibName} BEFORE PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${HOMME_Fortran_FLAGS_LIST}>)

    TARGET_INCLUDE_DIRECTORIES(
        ${dynLibName} PUBLIC
        ${SCREAM_INCLUDE_DIRS}
        ${EXEC_INCLUDE_DIRS}        # This comes from PREQX_KOKKOS_SETUP
        ${dynLibDir}                # For config.h
    )
    TARGET_INCLUDE_DIRECTORIES(
        ${dynLibName} SYSTEM PUBLIC
        ${SCREAM_TPL_INCLUDE_DIRS}
    )
    SET_TARGET_PROPERTIES (${dynLibName}
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${dynLibDir}"
      LIBRARY_OUTPUT_DIRECTORY "${dynLibDir}"
      RUNTIME_OUTPUT_DIRECTORY "${dynLibDir}"
    )
  ENDIF()
ENDMACRO(CreateDynamicsLib)

ADD_SUBDIRECTORY (tests)
