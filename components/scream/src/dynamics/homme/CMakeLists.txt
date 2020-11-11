# Set cmake config options for Homme
set (HOMME_SOURCE_DIR ${SCREAM_SOURCE_DIR}/../homme CACHE INTERNAL "Homme source directory")
set (HOMME_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/homme CACHE INTERNAL "Homme binary directory")

# If using Intel, we need to tell Homme to link against mkl rather than lapack

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  option (HOMME_USE_MKL "Whether to use Intel's MKL instead of blas/lapack" ON)
  option (HOMME_FIND_BLASLAPACK "Whether to use system blas/lapack" OFF)
else ()
  option (HOMME_USE_MKL "Whether to use Intel's MKL instead of blas/lapack" OFF)
  option (HOMME_FIND_BLASLAPACK "Whether to use system blas/lapack" ON)
endif ()

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
set(BUILD_HOMME_TOOL         OFF CACHE BOOL "")
# We DON'T want homme's pio support, so force this option
set(BUILD_HOMME_WITHOUT_PIOLIBRARY ON  CACHE BOOL "" FORCE)

# Only enable openMP in homme if Kokkos has openmp
option (ENABLE_OPENMP "OpenMP across elements" ${Kokkos_ENABLE_OPENMP})

# Enable targets in HOMME. For SCREAM, just enable theta, for now
set (BUILD_HOMME_THETA_KOKKOS ON CACHE BOOL "")

# This prevents Homme from building kokkos internally
# Note: if/when Homme will rely on EKAT for kokkos, you can get rid of this
set (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)

# configure Homme folder
add_subdirectory(${HOMME_SOURCE_DIR} ${HOMME_BINARY_DIR})

# Store homme compiler/linker flags, for future use, to prevent them from  going out of scope
# Long story: Homme sets all flags in CMAKE_<LANG>_FLAGS. Those are local vars, so at the end
#             of its main CMakeLists.txt, we copy them in the parent scope. But that's still
#             a local var, and the macro below will be likely called from a different subdir.
#             Therefore, store the flags internally in the cache, for future use in different
#             parts of the src tree.
# Also, strip leading/trailing spaces

string (STRIP "${HOMME_C_FLAGS}"       HOMME_C_FLAGS       )
string (STRIP "${HOMME_CXX_FLAGS}"     HOMME_CXX_FLAGS     )
string (STRIP "${HOMME_Fortran_FLAGS}" HOMME_Fortran_FLAGS )
string (STRIP "${HOMME_LINKER_FLAGS}"  HOMME_LINKER_FLAGS  )
set (HOMME_C_FLAGS ${HOMME_C_FLAGS} CACHE INTERNAL "")
set (HOMME_CXX_FLAGS ${HOMME_CXX_FLAGS} CACHE INTERNAL "")
set (HOMME_Fortran_FLAGS ${HOMME_Fortran_FLAGS} CACHE INTERNAL "")
set (HOMME_LINKER_FLAGS ${HOMME_LINKER_FLAGS} CACHE INTERNAL "")

#######################################################################
#     Macro to create dyn lib with specific compile-time options      #
#######################################################################

# This internal var allows to keep track of already created libs.
# If a test needs a lib already created, we simply retrieve it,
# otherwise we build a new library.
set (DynamicsLibsCreated   "" CACHE INTERNAL "List of created dynamics libraries names")

# NOTE: this macro creates TWO libraries: a homme library,
#       and a scream dynamics one.
macro (CreateDynamicsLib HOMME_TARGET NP PLEV QSIZE)
  if (NOT ${HOMME_TARGET} STREQUAL "theta-l_kokkos")
    message ("HOMME TARGET: ${HOMME_TARGET}")
    message (FATAL_ERROR "Error! So far, SCREAM only supports 'theta-l_kokkos' as dynamics target.")
  endif()

  # Create the lib name, and proceed building only if not already built
  set (hommeLibName "${HOMME_TARGET}_${NP}_${PLEV}_${QSIZE}")
  if ("${hommeLibName}" IN_LIST DynamicsLibsCreated)
    # This dynamics lib was built already somewhere in the project. Nothing to do
    set (dynLibName scream_${hommeLibName})
  else ()

    # Add this library name to the list of created libraries
    # Note: list(APPEND var items) does *not* work for cache variables
    set (DynamicsLibsCreated
         ${DynamicsLibsCreated}
         ${hommeLibName}
         CACHE INTERNAL "List of created dynamics libraries names"
    )
    # Set the dynamics library directory
    set (dynLibDir ${SCREAM_BINARY_DIR}/src/dynamics/homme/${hommeLibName})

    set (PREQX_NP ${NP})
    set (PREQX_PLEV ${PLEV})
    set (PREQX_QSIZE_D ${QSIZE})
    set (PREQX_USE_ENERGY FALSE)

    #####################################
    ###      Build Homme library      ###
    #####################################

    # Gather sources
    THETAL_KOKKOS_SETUP()

    # Re-set this, so that it is visible from wherever the macro is called
    set (HOMME_USE_KOKKOS TRUE)

    # Create the library
    createExecLib (${hommeLibName} ${HOMME_TARGET} "${THETAL_DEPS}" "${EXEC_LIB_INCLUDE_DIRS}"
                   ${PREQX_NP} ${PREQX_PLEV} ${PREQX_USE_ENERGY} ${QSIZE_D})
    target_compile_definitions(${hommeLibName} PUBLIC HOMMEXX_CONFIG_IS_CMAKE SCREAM)
    set_target_properties (${hommeLibName}
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${dynLibDir}"
      LIBRARY_OUTPUT_DIRECTORY "${dynLibDir}"
      RUNTIME_OUTPUT_DIRECTORY "${dynLibDir}"
    )

    if(${CMAKE_VERSION} VERSION_GREATER "3.13.0")
      target_link_options(${hommeLibName} PRIVATE "${HOMME_LINKER_FLAGS}")
    else ()
      target_link_libraries(${hommeLibName} "${HOMME_LINKER_FLAGS}")
    endif ()
    target_compile_options (${hommeLibName} PRIVATE $<$<COMPILE_LANGUAGE:C>:${HOMME_C_FLAGS}>)
    target_compile_options (${hommeLibName} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${HOMME_CXX_FLAGS}>)
    target_compile_options (${hommeLibName} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${HOMME_Fortran_FLAGS}>)

    #####################################
    ###     Build SCREAM library      ###
    #####################################

    # Gather sources
    set (SCREAM_DYNAMICS_SRC_DIR ${SCREAM_SRC_DIR}/dynamics/homme)

    set (SCREAM_DYNAMICS_SOURCES
      ${SCREAM_DYNAMICS_SRC_DIR}/atmosphere_dynamics.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/dynamics_driven_grids_manager.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_context_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_driver_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_grid_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_params_mod.F90
    )
    if (NOT SCREAM_LIB_ONLY)
      # We will also have unit tests, so add inputs initialization stuff
      list (APPEND SCREAM_DYNAMICS_SOURCES
        ${SCREAM_DYNAMICS_SRC_DIR}/homme_inputs_initializer.cpp
        ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_inputs_mod.F90
      )
    endif ()

    # Create library
    set (dynLibName scream_${hommeLibName})
    add_library(${dynLibName} ${SCREAM_DYNAMICS_SOURCES})
    target_link_libraries(${dynLibName} scream_share ${hommeLibName})

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
