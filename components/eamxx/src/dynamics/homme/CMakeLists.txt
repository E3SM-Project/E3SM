include (EkatSetCompilerFlags)

# Set cmake config options for Homme
set (HOMME_SOURCE_DIR ${SCREAM_BASE_DIR}/../homme CACHE INTERNAL "Homme source directory")
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
set(HOMME_ENABLE_COMPOSE     ON  CACHE BOOL "")
set(BUILD_HOMME_TOOL         OFF CACHE BOOL "")

if (NOT Kokkos_ENABLE_SERIAL)
  # We don't really *need* composef90, but Homme always builds it, so we need
  # to make sure it can build without errors.
  message(FATAL_ERROR "Homme's composef90 library requires Kokkos_ENABLE_SERIAL=ON.")
endif()

# We need homme's pio support for RRM grids, so force this option if this is a
# CIME build.
set(PIO_OFF TRUE)
if (SCREAM_CIME_BUILD)
  set(PIO_OFF FALSE)
endif ()
set(BUILD_HOMME_WITHOUT_PIOLIBRARY ${PIO_OFF} CACHE BOOL "" FORCE)

# Only enable openMP in homme if Kokkos has openmp
option (ENABLE_OPENMP "OpenMP across elements" ${Kokkos_ENABLE_OPENMP})

# Enable targets in HOMME. For SCREAM, just enable theta, for now
set (BUILD_HOMME_THETA_KOKKOS ON CACHE BOOL "")

# This prevents Homme from building kokkos internally
# Note: if/when Homme will rely on EKAT for kokkos, you can get rid of this
set (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)

# Make Homme use the same pack size that scream does
set (HOMMEXX_VECTOR_SIZE ${SCREAM_PACK_SIZE} CACHE STRING "")

# Set Homme mpi on device option
option (HOMMEXX_MPI_ON_DEVICE "Whether we want to use device pointers for MPI calls (relevant only for GPU builds)" ${SCREAM_MPI_ON_DEVICE})

# configure Homme folder
add_subdirectory(${HOMME_SOURCE_DIR} ${HOMME_BINARY_DIR})

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
    message (FATAL_ERROR "Error! So far, SCREAM only supports 'theta-l_kokkos' as dynamics target.")
  endif()

  add_definitions(-DHAVE_CONFIG_H)

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
    # Link to cime's csm_share lib
    target_link_libraries (${hommeLibName} csm_share)
    target_link_libraries (${hommeLibName} composec++)
    SetCudaFlags(${hommeLibName})
    SetOmpFlags(${hommeLibName})

    #####################################
    ###     Build SCREAM library      ###
    #####################################

    # Gather sources
    set (SCREAM_DYNAMICS_SRC_DIR ${SCREAM_SRC_DIR}/dynamics/homme)

    set (SCREAM_DYNAMICS_SOURCES
      ${SCREAM_DYNAMICS_SRC_DIR}/atmosphere_dynamics.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/atmosphere_dynamics_fv_phys.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/atmosphere_dynamics_rayleigh_friction.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/physics_dynamics_remapper.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/homme_grids_manager.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_context_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_driver_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_grid_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/homme_params_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/dyn_grid_mod.F90
      ${SCREAM_DYNAMICS_SRC_DIR}/interface/phys_grid_mod.F90
    )

    # Create library
    set (dynLibName scream_${hommeLibName})
    add_library(${dynLibName} ${SCREAM_DYNAMICS_SOURCES})
    target_link_libraries(${dynLibName} scream_share scream_io ${hommeLibName})
    get_target_property(modulesDir ${hommeLibName} Fortran_MODULE_DIRECTORY)
    set_target_properties(${dynLibName} PROPERTIES Fortran_MODULE_DIRECTORY ${modulesDir})
    target_include_directories(${dynLibName} PUBLIC ${modulesDir})
  endif()
endmacro(CreateDynamicsLib)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory (tests)
endif()
