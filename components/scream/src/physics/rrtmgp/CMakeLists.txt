include(EkatUtils)
include(EkatSetCompilerFlags)
include(ScreamUtils)

# Add CUDA flags for YAKL
set(YAKL_CXX_FLAGS "")
if (CUDA_BUILD)
  set(YAKL_ARCH "CUDA")
  list(APPEND YAKL_CXX_FLAGS -DYAKL_ARCH_CUDA --expt-extended-lambda --expt-relaxed-constexpr)
endif()

# RRTMGP++ requires YAKL
add_subdirectory(${SCREAM_BASE_DIR}/../../externals/YAKL ${CMAKE_BINARY_DIR}/externals/YAKL)
target_compile_options(yakl PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${YAKL_CXX_FLAGS}>)
EkatDisableAllWarning(yakl)

# For debug builds, set -DYAKL_DEBUG
string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_ci)
if (CMAKE_BUILD_TYPE_ci STREQUAL "debug")
  target_compile_definitions(yakl PUBLIC YAKL_DEBUG)
endif()

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
SetCudaFlags(rrtmgp)
target_link_libraries(rrtmgp PUBLIC yakl)
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
set(INTERFACE_SRC
  atmosphere_radiation.cpp
  scream_rrtmgp_interface.cpp
  mo_garand_atmos_io.cpp
  mo_load_coefficients.cpp
  mo_load_cloud_coefficients.cpp
  rrtmgp_test_utils.cpp
  shr_orb_mod_c2f.F90
)
add_library(scream_rrtmgp ${INTERFACE_SRC})
find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATH}/lib)
target_link_libraries(scream_rrtmgp PUBLIC ${NETCDF_C} rrtmgp yakl scream_share physics_share csm_share)
target_include_directories(scream_rrtmgp SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${NetCDF_C_PATH}/include
    ${EAM_RRTMGP_DIR}/external)

# Ensure RRTMGP lookup tables are present in the data dir
set (RRTMGP_TABLES
  init/rrtmgp-data-sw-g224-2018-12-04.nc
  init/rrtmgp-data-lw-g256-2018-12-04.nc
  init/rrtmgp-allsky.nc
  init/rrtmgp-cloud-optics-coeffs-sw.nc
  init/rrtmgp-cloud-optics-coeffs-lw.nc
)

foreach (file IN ITEMS ${RRTMGP_TABLES})
  GetInputFile(${file})
endforeach()

# Build tests
if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
