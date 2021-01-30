# YAKL requires netcdf-cxx4
#execute_process(COMMAND nc-config --libdir OUTPUT_VARIABLE NETCDF_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
#message(STATUS "NETCDF_LIBDIR: ${NETCDF_LIBDIR}")
#set(netCDF_LIBRARIES ${NETCDF_LIBDIR}/libnetcdf.so)
#execute_process(COMMAND nc-config --includedir OUTPUT_VARIABLE netCDF_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
#message(STATUS "netCDF_LIBRARIES: ${netCDF_LIBRARIES}")
#message(STATUS "netCDF_INCLUDE_DIR: ${netCDF_INCLUDE_DIR}")
#add_subdirectory(${SCREAM_BASE_DIR}/../../externals/netcdf-cxx4 ${CMAKE_CURRENT_BINARY_DIR}/netcdf-cxx4_build)
#target_include_directories(yakl PUBLIC ${SCREAM_BASE_DIR}/../../externals/netcdf-cxx4/cxx4)

# RRTMGP++ requires YAKL
add_subdirectory(${SCREAM_BASE_DIR}/../../externals/YAKL ${CMAKE_CURRENT_BINARY_DIR}/YAKL_build)
target_compile_options(yakl PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fno-default-real-8 -fno-default-double-8>)

# Add RRTMGP source files.
set(RRTMGP_SRCS
  atmosphere_radiation.cpp
  scream_rrtmgp_interface.cpp
  rrtmgp_inputs_initializer.cpp
  mo_garand_atmos_io.cpp
  mo_load_coefficients.cpp
  mo_load_cloud_coefficients.cpp
  tests/rrtmgp_test_utils.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels/mo_gas_optics_kernels.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/mo_rrtmgp_constants.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/mo_rrtmgp_util_reorder.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/mo_rte_lw.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/mo_rte_sw.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/expand_and_transpose.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels/mo_fluxes_broadband_kernels.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels/mo_optical_props_kernels.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels/mo_rte_solver_kernels.cpp
)

set(RRTMGP_HEADERS
  rrtmgp.hpp
  scream_rrtmgp_interface.hpp
  atmosphere_radiation.hpp
  rrtmgp_inputs_initializer.hpp
)

add_library(rrtmgp ${RRTMGP_SRCS})
set_source_files_properties(${RRTMGP_SRCS} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
target_include_directories(rrtmgp PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(rrtmgp PUBLIC ./)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/examples)
# Again, this is a hack to get the standalone tests to build that need the mo_garand_atmos_io.h header
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/examples/all-sky)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/extensions/cloud_optics)
target_include_directories(rrtmgp PUBLIC ${SCREAM_BASE_DIR}/../../externals/YAKL)
target_include_directories(rrtmgp PUBLIC tests)
set_target_properties(rrtmgp PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
set_target_properties(rrtmgp PROPERTIES CXX_STANDARD 14)
target_link_libraries(rrtmgp yakl physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
