# Set cmake config options for Homme
set (HOMME_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/homme CACHE INTERNAL "Homme source directory")
set (HOMME_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/homme CACHE INTERNAL "Homme binary directory")

# Set some important options used to configure homme
set (WITH_PNETCDF FALSE CACHE BOOL "Whether pnetcdf is available.")
set (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "Whether homme should find (rther than build) blas/lapack")

# Disable all the targets by default
set(BUILD_HOMME_SWEQX        OFF CACHE BOOL "")
set(BUILD_HOMME_PREQX        OFF CACHE BOOL "")
set(BUILD_HOMME_THETA        OFF CACHE BOOL "")
set(BUILD_HOMME_PREQX_ACC    OFF CACHE BOOL "")
set(BUILD_HOMME_PREQX_KOKKOS OFF CACHE BOOL "")
set(BUILD_HOMME_PESE         OFF CACHE BOOL "")
set(BUILD_HOMME_SWIM         OFF CACHE BOOL "")
set(BUILD_HOMME_PRIM         OFF CACHE BOOL "")
set(HOMME_ENABLE_COMPOSE     OFF CACHE BOOL "")
set(HOMME_FIND_BLASLAPACK    ON  CACHE BOOL "")
set(BUILD_HOMME_WITHOUT_PIOLIBRARY ON  CACHE BOOL "")
set(BUILD_HOMME_TOOL         OFF CACHE BOOL "")

if (CIME_BUILD)
  # Test vars from Macros.cmake
  if (NETCDF_C_PATH)
    if (NOT NETCDF_FORTRAN_PATH)
      message(FATAL_ERROR "NETCDF_C_PATH specified without NETCDF_FORTRAN_PATH")
    endif()
    set(NetCDF_Fortran_PATH ${NETCDF_FORTRAN_PATH} CACHE PATH "Netcdf fortran lib path")
  elseif()
    if (NETCDF_FORTRAN_PATH)
      message(FATAL_ERROR "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")
    endif()
  else()
    set(NetCDF_Fortran_PATH ${NETCDF_PATH} CACHE PATH "Netcdf fortran lib path")
  endif()
endif()

# List of acceptable homme targets for scream
# Only these targets will be usable from scream
set (HOMME_TARGETS
     preqx_kokkos
     # THETA_KOKKOS # coming soon...
)

# Enable the requested target
foreach (target IN ITEMS ${HOMME_TARGETS})
  if ("${target}" STREQUAL "preqx_kokkos")
    set (BUILD_HOMME_PREQX_KOKKOS TRUE)
  elseif ("${target}" STREQUAL "theta_kokkos")
    set (BUILD_HOMME_THETA_KOKKOS TRUE)
  else ()
    messagE (FATAL_ERROR "Homme target '${target}' is not supported within scream.\n")
  endif()
endforeach()

set (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)
add_subdirectory(homme)

# Take the variable that homme set, and mark them as (internal) cached, so they are globally visible.
# This is needed since the macro CreateDynamicsLib will be invoked from outside this directory,
# for instance, from the scream/tests folder.
SET (HOMME_C_FLAGS       ${HOMME_C_FLAGS}       CACHE INTERNAL "List of C compiler flags for Homme")
SET (HOMME_CXX_FLAGS     ${HOMME_CXX_FLAGS}     CACHE INTERNAL "List of CXX compiler flags for Homme")
SET (HOMME_Fortran_FLAGS ${HOMME_Fortran_FLAGS} CACHE INTERNAL "List of Fortran compiler flags for Homme")

SET (DynamicsLibsCreated   "" CACHE INTERNAL "List of created dynamics libraries names")

macro (CreateDynamicsLib HOMME_TARGET NP PLEV QSIZE USE_PIO USE_ENERGY)
  set (PREQX_NP ${NP})
  set (PREQX_PLEV ${PLEV})
  set (PREQX_QSIZE_D ${QSIZE})
  set (PREQX_USE_PIO ${USE_PIO})
  set (PREQX_USE_ENERGY ${USE_ENERGY})

  set (dynLibName "${HOMME_TARGET}_np${NP}_nlev${PLEV}_qsize${QSIZE}_pio${USE_PIO}_energy${USE_ENERGY}")
  if ("${dynLibName}" IN_LIST DynamicsLibsCreated)
    set (dynLibDir ${SCREAM_BINARY_DIR}/src/dynamics/homme/${dynLibName})
  else ()
    # Add this library name to the list of created libraries
    set (DynamicsLibsCreated
         ${DynamicsLibsCreated}
         ${dynLibName}
         CACHE INTERNAL "List created dynamics libraries names"
    )

    # Set the dynamics library directory
    set (dynLibDir ${SCREAM_BINARY_DIR}/src/dynamics/homme/${dynLibName})

    add_definitions(-DHOMMEXX_CONFIG_IS_CMAKE)
    set (SCREAM_DYNAMICS_SRC_DIR ${SCREAM_SRC_DIR}/dynamics/homme)

    # re-set these, so that they are visible from wherever the macro is called
    set (HOMME_SOURCE_DIR ${SCREAM_SOURCE_DIR}/src/dynamics/homme/homme)
    set (HOMME_BINARY_DIR ${SCREAM_BINARY_DIR}/src/dynamics/homme/homme)

    if ("${HOMME_TARGET}" STREQUAL "preqx_kokkos")
      # Setup the homme target
      PREQX_KOKKOS_SETUP()
    else ()
      message (FATAL_ERROR "Homme target '${HOMME_TARGET}' is not supported within scream.\n \
                            Did you forget to set 'HOMME_TARGET' before calling 'CreateDynamicsLib'?\n")
    endif()

    set (SCREAM_DYNAMICS_SOURCES
         ${SCREAM_DYNAMICS_SRC_DIR}/atmosphere_dynamics.cpp
         ${SCREAM_DYNAMICS_SRC_DIR}/dynamics_driven_grids_manager.cpp
         ${SCREAM_DYNAMICS_SRC_DIR}/scream_homme_interface.F90
         ${PREQX_DEPS_CXX}
         ${PREQX_DEPS_F90}
    )

    # append all libraries need byhomme to the list of libraries needed by dynamics
    set (SCREAM_DYNAMICS_LINK_LIBRARIES
         scream_share
         timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
    )
    if (NOT BUILD_HOMME_WITHOUT_PIOLIBRARY)
      message("************adding pio")
      set (SCREAM_DYNAMICS_LINK_LIBRARIES
           ${SCREAM_DYNAMICS_LINK_LIBRARIES}
           pio
      )
    endif()

    # Set the variables that need to be expanded by cmake during configure_file
    set (NUM_POINTS ${PREQX_NP})
    set (NUM_CELLS ${PREQX_NC})
    set (NUM_PLEV ${PREQX_PLEV})
    set (QSIZE_D ${PREQX_QSIZE_D})

    set (PIO ${PREQX_USE_PIO})
    set (PIO_INTERP NOT ${PREQX_USE_PIO})
    set (ENERGY_DIAGNOSTICS ${PREQX_USE_ENERGY})

    # this is needed to compile the test executables with the correct options for homme
    set(THIS_CONFIG_IN ${HOMME_SOURCE_DIR}/src/${HOMME_TARGET}/config.h.cmake.in)
    set(THIS_CONFIG_HC ${dynLibDir}/config.h.c)
    set(THIS_CONFIG_H ${dynLibDir}/config.h)

    # first configure the file (which formats the file as C)
    hommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

    add_definitions(-DHAVE_CONFIG_H)

    # convert strings to lists
    # note: TARGET_COMPILE_OPTIONS perform 'de-duplication' of repeated flags.
    #       this causes the repeated "-fp-model X -fp-model Y" to be turned into
    #       "-fp-model X Y", which the compiler does not like.
    #       Thankfully, there is a remedy, albeit a bit ugly: prepend the flag to be kept
    #       (even if repeated) with 'SHELL: '.
    string(REPLACE "-fp-model" "SHELL: -fp-model" HOMME_C_FLAGS_LIST_TMP       "${HOMME_C_FLAGS}")
    string(REPLACE "-fp-model" "SHELL: -fp-model" HOMME_CXX_FLAGS_LIST_TMP     "${HOMME_CXX_FLAGS}")
    string(REPLACE "-fp-model" "SHELL: -fp-model" HOMME_Fortran_FLAGS_LIST_TMP "${HOMME_Fortran_FLAGS}")
    string(REPLACE " " ";" HOMME_C_FLAGS_LIST       "${HOMME_C_FLAGS_TMP}")
    string(REPLACE " " ";" HOMME_CXX_FLAGS_LIST     "${HOMME_CXX_FLAGS_TMP}")
    string(REPLACE " " ";" HOMME_Fortran_FLAGS_LIST "${HOMME_Fortran_FLAGS_TMP}")

    add_library(${dynLibName} ${SCREAM_DYNAMICS_SOURCES})
    target_link_libraries(${dynLibName} ${SCREAM_DYNAMICS_LINK_LIBRARIES})
    set_target_properties (${dynLibName} PROPERTIES Fortran_MODULE_DIRECTORY ${dynLibDir}/${dynLibName}_modules)

    # Add flags, using a generator to achieve language specific behavior
    target_compile_options(${dynLibName} BEFORE PUBLIC $<$<COMPILE_LANGUAGE:C>:${HOMME_C_FLAGS_LIST}>)
    target_compile_options(${dynLibName} BEFORE PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${HOMME_CXX_FLAGS_LIST}>)
    target_compile_options(${dynLibName} BEFORE PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${HOMME_Fortran_FLAGS_LIST}>)

    target_include_directories(
        ${dynLibName} PUBLIC
        ${SCREAM_INCLUDE_DIRS}
        ${EXEC_INCLUDE_DIRS}        # This comes from PREQX_KOKKOS_SETUP
        ${dynLibDir}                # For config.h
    )
    target_include_directories(
        ${dynLibName} SYSTEM PUBLIC
        ${SCREAM_TPL_INCLUDE_DIRS}
    )
    set_target_properties (${dynLibName}
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${dynLibDir}"
      LIBRARY_OUTPUT_DIRECTORY "${dynLibDir}"
      RUNTIME_OUTPUT_DIRECTORY "${dynLibDir}"
    )
  endif()
endmacro(CreateDynamicsLib)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory (tests)
endif()
