# Add CUDA flags for YAKL
if (CUDA_BUILD)
    set(ARCH "CUDA")
    set(YAKL_CXX_FLAGS "-D__USE_CUDA__ --expt-extended-lambda --expt-relaxed-constexpr ${CUDA_FLAGS} ${YAKL_CXX_FLAGS}")
endif()

# RRTMGP++ requires YAKL
add_subdirectory(${SCREAM_BASE_DIR}/../../externals/YAKL ${CMAKE_CURRENT_BINARY_DIR}/YAKL_build)
target_compile_options(yakl PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fno-default-real-8 -fno-default-double-8>)

# Build RRTMGP library; this builds the core RRTMGP external source as a library named "rrtmgp"
# NOTE: The external RRTMGP build needs some fixes to work with CUDA in a library build, so for now we will build these ourselves
#add_subdirectory(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp ${CMAKE_CURRENT_BINARY_DIR}/rrtmgp_build)
set(EXTERNAL_SRC
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels/mo_gas_optics_kernels.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/mo_rrtmgp_constants.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/mo_rrtmgp_util_reorder.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/expand_and_transpose.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels/mo_fluxes_broadband_kernels.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels/mo_optical_props_kernels.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels/mo_rte_solver_kernels.cpp
)
add_library(rrtmgp ${EXTERNAL_SRC})
target_link_libraries(rrtmgp yakl)
set_target_properties(rrtmgp PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS} ${CXX_FLAGS}")
target_include_directories(rrtmgp PUBLIC 
    ${SCREAM_BASE_DIR}/../../externals/YAKL
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels
)

# Build RRTMGP interface; note that we separate the SCREAM-specific RRTMGP interface
# from the external core RRTMGP library because, ideally, the RRTMGP library has its
# own build, and we would just use add_subdirectory() above to build it. Also, this
# separates out the code that comprises the core RRTMGP library from the extensions
# and examples that we have modified for use in SCREAM specifically.
set(INTERFACE_SRC
  atmosphere_radiation.cpp
  scream_rrtmgp_interface.cpp
  rrtmgp_inputs_initializer.cpp
  mo_garand_atmos_io.cpp
  mo_load_coefficients.cpp
  mo_load_cloud_coefficients.cpp
  tests/rrtmgp_test_utils.cpp
)
add_library(scream_rrtmgp ${INTERFACE_SRC})
set_target_properties(scream_rrtmgp PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
target_link_libraries(scream_rrtmgp ${NETCDF_C} rrtmgp yakl scream_share physics_share)
target_include_directories(scream_rrtmgp PUBLIC ${NetCDF_C_PATHS}/include)
target_include_directories(scream_rrtmgp PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(scream_rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../../externals)
target_include_directories(scream_rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external)

# Build tests
if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
