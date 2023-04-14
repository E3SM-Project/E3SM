include(EkatUtils)
include(EkatSetCompilerFlags)
include(ScreamUtils)

##################################
#             YAKL               #
##################################

# RRTMGP++ requires YAKL

if (TARGET yakl)
  # Other E3SM components are building YAKL...
  message ("It appears some other part of E3SM is building YAKL.\n"
           "We will reuse that, but if this is a debug build we will\n"
           "add the --fmad=false flag to the cuda flags used by YAKL\n")
else ()
  # Prepare CUDA/HIP flags for YAKL
  if (CUDA_BUILD)
    set(YAKL_ARCH "CUDA")
    set(YAKL_CUDA_FLAGS "-DYAKL_ARCH_CUDA --expt-extended-lambda --expt-relaxed-constexpr -ccbin ${CMAKE_CXX_COMPILER}")
    string (REPLACE " " ";" YAKL_CUDA_FLAGS_LIST ${YAKL_CUDA_FLAGS})
  endif()
  if (HIP_BUILD)
    set(YAKL_ARCH "HIP")
    set(YAKL_HIP_FLAGS "-DYAKL_ARCH_HIP -O3 -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a -x hip")
    string (REPLACE " " ";" YAKL_HIP_FLAGS_LIST ${YAKL_HIP_FLAGS})
  endif()

  set (YAKL_SOURCE_DIR ${SCREAM_BASE_DIR}/../../externals/YAKL)
  add_subdirectory(${YAKL_SOURCE_DIR} ${CMAKE_BINARY_DIR}/externals/YAKL)

  # Set some additional flag/cpp option on the yakl target

  cmake_policy (SET CMP0079 NEW) # Allow to link to a tgt from a different directory

  # EAMxx *requires* MPI, so simply look for it, then link against it
  find_package(MPI REQUIRED COMPONENTS C)
  target_link_libraries (yakl MPI::MPI_C)
  EkatDisableAllWarning(yakl)

  # For debug builds, set -DYAKL_DEBUG
  string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_ci)
  if (CMAKE_BUILD_TYPE_ci STREQUAL "debug")
    target_compile_definitions(yakl PUBLIC YAKL_DEBUG)
  endif()
endif()

# Whether we built yakl or it was already built, we still need to ensure
# debug builds have --fmad=false in it. The EKAT util automatically sets
# the flag for debug builds
SetCudaFlags(yakl CUDA_LANG)

##################################
#           RRTMGP               #
##################################

list(APPEND CMAKE_MODULE_PATH ${YAKL_SOURCE_DIR})
include (yakl_utils)

set(EAM_RRTMGP_DIR ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp)
# Build RRTMGP library; this builds the core RRTMGP external source as a library named "rrtmgp"
# NOTE: The external RRTMGP build needs some fixes to work with CUDA in a library build, so for now we will build these ourselves
set(EXTERNAL_SRC
  ${EAM_RRTMGP_DIR}/external/cpp/rrtmgp/kernels/mo_gas_optics_kernels.cpp
  ${EAM_RRTMGP_DIR}/external/cpp/rrtmgp/mo_rrtmgp_constants.cpp
  ${EAM_RRTMGP_DIR}/external/cpp/rrtmgp/mo_rrtmgp_util_reorder.cpp
  ${EAM_RRTMGP_DIR}/external/cpp/rte/expand_and_transpose.cpp
  ${EAM_RRTMGP_DIR}/external/cpp/rte/kernels/mo_fluxes_broadband_kernels.cpp
  ${EAM_RRTMGP_DIR}/external/cpp/rte/kernels/mo_optical_props_kernels.cpp
  ${EAM_RRTMGP_DIR}/external/cpp/rte/kernels/mo_rte_solver_kernels.cpp
  ${EAM_RRTMGP_DIR}/external/cpp/extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
)
add_library(rrtmgp ${EXTERNAL_SRC})
EkatDisableAllWarning(rrtmgp)
yakl_process_target(rrtmgp)

# NOTE: cannot use 'PUBLIC' in target_link_libraries,
#       since yakl_process_target already used it
#       with the "plain" signature
target_link_libraries(rrtmgp yakl)
target_include_directories(rrtmgp PUBLIC
    ${SCREAM_BASE_DIR}/../../externals/YAKL
    ${EAM_RRTMGP_DIR}/external/cpp
    ${EAM_RRTMGP_DIR}/external/cpp/rte
    ${EAM_RRTMGP_DIR}/external/cpp/rte/kernels
    ${EAM_RRTMGP_DIR}/external/cpp/rrtmgp
    ${EAM_RRTMGP_DIR}/external/cpp/rrtmgp/kernels
)

# The lines below are needed to ensure that kokkos_launch_compiler injects
# nvcc into compilations. rrtmgp uses YAKL, not kokkos, so the wrapper
# didn't know to add nvcc without these lines.
target_compile_definitions(rrtmgp PRIVATE KOKKOS_DEPENDENCE)
target_link_options(rrtmgp PRIVATE -DKOKKOS_DEPENDENCE)

# Build RRTMGP interface; note that we separate the SCREAM-specific RRTMGP interface
# from the external core RRTMGP library because, ideally, the RRTMGP library has its
# own build, and we would just use add_subdirectory() above to build it. Also, this
# separates out the code that comprises the core RRTMGP library from the extensions
# and examples that we have modified for use in SCREAM specifically.

# However, due to the mix of YAKL and Kokkos, we split the target in two:
#  - scream_rrtmgp: kokkos-based interface to EAMxx
#  - scream_rrtmgp_yakl: source codes to be built with YAKL flags/options

##################################
#       SCREAM_RRTMGP_YAKL       #
##################################

set(SCREAM_RRTMGP_SOURCES_YAKL
  scream_rrtmgp_interface.cpp
  mo_garand_atmos_io.cpp
  mo_load_coefficients.cpp
  mo_load_cloud_coefficients.cpp
)

add_library(scream_rrtmgp_yakl ${SCREAM_RRTMGP_SOURCES_YAKL})
yakl_process_target(scream_rrtmgp_yakl)

# NOTE: cannot use 'PUBLIC' in target_link_libraries,
#       since yakl_process_target already used it
#       with the "plain" signature
find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATH}/lib)
target_link_libraries(scream_rrtmgp_yakl ${NETCDF_C} rrtmgp scream_share)
target_include_directories(scream_rrtmgp_yakl PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(scream_rrtmgp_yakl SYSTEM PUBLIC
    ${NetCDF_C_PATH}/include
    ${EAM_RRTMGP_DIR}/external)

##################################
#        SCREAM_RRTMGP           #
##################################

set(SCREAM_RRTMGP_SOURCES
  atmosphere_radiation.cpp
  shr_orb_mod_c2f.F90
)

add_library(scream_rrtmgp ${SCREAM_RRTMGP_SOURCES})
target_link_libraries(scream_rrtmgp PUBLIC scream_share physics_share csm_share scream_rrtmgp_yakl)
set_target_properties(scream_rrtmgp PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(scream_rrtmgp PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/modules)

# If yakl builds with LANG!=CXX, then the yakl CPP defines don't transfer to scream
# targets, b/c of the lang difference. So, if YAKL_ARCH is set, we add
# ${YAKL_${YAKL_ARCH}_FLAGS} flags to the CXX flags of scream_rrtmgp.
# In particular, this will ensure that all the yakl macros
# are correctly defined in YAKL headers, depending on the backend
if (YAKL_ARCH)
  target_compile_options(scream_rrtmgp PUBLIC
        "$<$<COMPILE_LANGUAGE:CXX>:${YAKL_${YAKL_ARCH}_FLAGS_LIST}>")
endif()


# Ensure RRTMGP lookup tables are present in the data dir
set (RRTMGP_TABLES
  scream/init/rrtmgp-data-sw-g112-210809.nc
  scream/init/rrtmgp-data-lw-g128-210809.nc
  scream/init/rrtmgp-allsky.nc
  scream/init/rrtmgp-cloud-optics-coeffs-sw.nc
  scream/init/rrtmgp-cloud-optics-coeffs-lw.nc
)

foreach (file IN ITEMS ${RRTMGP_TABLES})
  GetInputFile(${file})
endforeach()

# Build tests
if (NOT SCREAM_LIB_ONLY)
  # Some utils for the tests. Build them only once, for all execs that might need them
  add_library (rrtmgp_test_utils rrtmgp_test_utils.cpp)
  target_link_libraries(rrtmgp_test_utils PUBLIC scream_rrtmgp)

  add_subdirectory(tests)
endif()
