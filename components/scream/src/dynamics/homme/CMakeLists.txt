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

# We DON'T want homme's pio support, so force this option
set(BUILD_HOMME_WITHOUT_PIOLIBRARY ON  CACHE BOOL "" FORCE)

# Only enable openMP in homme if Kokkos has openmp
option (ENABLE_OPENMP "OpenMP across elements" ${Kokkos_ENABLE_OPENMP})

# Enable targets in HOMME. For SCREAM, just enable theta, for now
set (BUILD_HOMME_THETA_KOKKOS ON CACHE BOOL "")

# This prevents Homme from building kokkos internally
# Note: if/when Homme will rely on EKAT for kokkos, you can get rid of this
set (E3SM_INTERNAL_KOKKOS_ALREADY_BUILT TRUE)

# Note: we need scream's fp model to superseed homme's one. So we append our
# current compiler flags to homme's flags
set (TEMP_C_FLAGS ${CMAKE_C_FLAGS})
set (TEMP_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set (TEMP_F_FLAGS ${CMAKE_Fortran_FLAGS})
set (CMAKE_C_FLAGS)
set (CMAKE_CXX_FLAGS)
set (CMAKE_Fortran_FLAGS)

# Make Homme use the same pack size that scream does
set (HOMMEXX_VECTOR_SIZE ${SCREAM_PACK_SIZE} CACHE STRING "")

# configure Homme folder
add_subdirectory(${HOMME_SOURCE_DIR} ${HOMME_BINARY_DIR})

# *APPEND* scream flags to homme flags
set (HOMME_C_FLAGS       ${HOMME_C_FLAGS} ${TEMP_C_FLAGS})
set (HOMME_CXX_FLAGS     ${HOMME_CXX_FLAGS} ${TEMP_CXX_FLAGS})
set (HOMME_Fortran_FLAGS ${HOMME_Fortran_FLAGS} ${TEMP_F_FLAGS})

# A bunch of string manipulation to get homme compiler/link flags correctly
# forwarded to our libs. Here are the takeaways:
#
# - we need to store HOMME_<LANG>_FLAGS in the cache, since CreateDynamicsLib
#   may be called from somewhere else in the src tree, and those vars would
#   otherwise get out of scope
# - We need to turn the strings into lists, by replacing ' ' with ';'
# - The hard one: target_compile_options performs 'de-duplication' of repeated flags,
#   meaning that '-D A -D B' is turned into '-D A B'. For fp model on Intel, this
#   would cause "-fp-model X -fp-model Y" to be turned into "-fp-model X Y",
#   which the compiler does not like. Thankfully, there is a remedy, albeit a bit ugly:
#   prepend the flag to be kept (even if repeated) with 'SHELL: '.

string (REPLACE " " ";" HOMME_C_FLAGS_TMP       "${HOMME_C_FLAGS}")
string (REPLACE " " ";" HOMME_CXX_FLAGS_TMP     "${HOMME_CXX_FLAGS}")
string (REPLACE " " ";" HOMME_Fortran_FLAGS_TMP "${HOMME_Fortran_FLAGS}")
string (REPLACE " " ";" HOMME_LINKER_FLAGS_TMP  "${HOMME_LINKER_FLAGS}")
string(REGEX REPLACE "-fp-model;(fast|fast=1|fast=2|strict|precise|source|consistent);" "SHELL:-fp-model \\1;" HOMME_C_FLAGS_TMP       "${HOMME_C_FLAGS_TMP}")
string(REGEX REPLACE "-fp-model;(fast|fast=1|fast=2|strict|precise|source|consistent);" "SHELL:-fp-model \\1;" HOMME_CXX_FLAGS_TMP     "${HOMME_CXX_FLAGS_TMP}")
string(REGEX REPLACE "-fp-model;(fast|fast=1|fast=2|strict|precise|source|consistent);" "SHELL:-fp-model \\1;" HOMME_Fortran_FLAGS_TMP "${HOMME_Fortran_FLAGS_TMP}")
set (HOMME_C_FLAGS "${HOMME_C_FLAGS_TMP}" CACHE INTERNAL "")
set (HOMME_CXX_FLAGS "${HOMME_CXX_FLAGS_TMP}" CACHE INTERNAL "")
set (HOMME_Fortran_FLAGS "${HOMME_Fortran_FLAGS_TMP}" CACHE INTERNAL "")
set (HOMME_LINKER_FLAGS "${HOMME_LINKER_FLAGS_TMP}" CACHE INTERNAL "")

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

    if (SCREAM_CIME_BUILD)
      ### HACK ALERT! ###

      # Remove shr_blah_mod.F90 from homme/utils/csm_share,
      # since we use those from $repo_root/share/util.
      # Note: the cleanest thing would be to mod homme, to use cime stuff
      set (REMOVE_FILES)
      foreach (file IN_LIST ${THETAL_DEPS})
        if (file MATCHES "csm_share/shr_")
          list (APPEND REMOVE_FILES ${file})
        endif ()
      endforeach()
      list (REMOVE_ITEM THETAL_DEPS ${REMOVE_FILES})
    endif ()

    # Re-set this, so that it is visible from wherever the macro is called
    set (HOMME_USE_KOKKOS TRUE)

    # Create the library
    createExecLib (${hommeLibName} ${HOMME_TARGET} "${THETAL_DEPS}" "${EXEC_LIB_INCLUDE_DIRS}"
                   ${PREQX_NP} ${PREQX_PLEV} ${PREQX_USE_ENERGY} ${QSIZE_D})
    target_compile_definitions(${hommeLibName} PUBLIC HOMMEXX_CONFIG_IS_CMAKE SCREAM)
    if (SCREAM_CIME_BUILD)
      # Link to cime's csm_share lib
      target_link_libraries (${hommeLibName} csm_share)
    endif ()
    target_link_libraries (${hommeLibName} composec++)
    SetCudaFlags(${hommeLibName})
    SetOmpFlags(${hommeLibName})

    # target_link_options is supported since 3.13, but there's a bug with static libs until 3.17
    # see https://gitlab.kitware.com/cmake/cmake/-/issues/20022 for details
    if(${CMAKE_VERSION} VERSION_GREATER "3.17.0")
      target_link_options(${hommeLibName} PUBLIC "${HOMME_LINKER_FLAGS}")
    else ()
      target_link_libraries(${hommeLibName} "${HOMME_LINKER_FLAGS}")
    endif ()
    target_compile_options (${hommeLibName} BEFORE PRIVATE $<$<COMPILE_LANGUAGE:C>:${HOMME_C_FLAGS}>)
    target_compile_options (${hommeLibName} BEFORE PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${HOMME_CXX_FLAGS}>)
    target_compile_options (${hommeLibName} BEFORE PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${HOMME_Fortran_FLAGS}>)

    #####################################
    ###     Build SCREAM library      ###
    #####################################

    # Gather sources
    set (SCREAM_DYNAMICS_SRC_DIR ${SCREAM_SRC_DIR}/dynamics/homme)

    set (SCREAM_DYNAMICS_SOURCES
      ${SCREAM_DYNAMICS_SRC_DIR}/atmosphere_dynamics.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/physics_dynamics_remapper.cpp
      ${SCREAM_DYNAMICS_SRC_DIR}/dynamics_driven_grids_manager.cpp
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
    target_link_libraries(${dynLibName} scream_share ${hommeLibName})
    get_target_property(modulesDir ${hommeLibName} Fortran_MODULE_DIRECTORY)
    set_target_properties(${dynLibName} PROPERTIES Fortran_MODULE_DIRECTORY ${modulesDir})
    target_include_directories(${dynLibName} PUBLIC ${modulesDir})
  endif()
endmacro(CreateDynamicsLib)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory (tests)
endif()
