include(EkatUtils)
include(EkatSetCompilerFlags)
include(ScreamUtils)

##################################
#           RRTMGP               #
##################################

set(EAM_RRTMGP_DIR ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp)
# Build RRTMGP library; this builds the core RRTMGP external source as a library named "rrtmgp"
# NOTE: The external RRTMGP build needs some fixes to work with CUDA in a library build, so for now we will build these ourselves
add_library(rrtmgp INTERFACE)
target_compile_definitions(rrtmgp INTERFACE EAMXX_HAS_RRTMGP)
if (Kokkos_ENABLE_CUDA)
  target_compile_options(rrtmgp INTERFACE $<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr>)
endif()

if (NOT TARGET Kokkos::kokkos)
  find_package(Kokkos REQUIRED)
endif ()
target_link_libraries(rrtmgp INTERFACE Kokkos::kokkos)
target_include_directories(rrtmgp INTERFACE
    ${EAM_RRTMGP_DIR}/external/cpp
    ${EAM_RRTMGP_DIR}/external/cpp/extensions/cloud_optics
    ${EAM_RRTMGP_DIR}/external/cpp/examples
    ${EAM_RRTMGP_DIR}/external/cpp/examples/all-sky
    ${EAM_RRTMGP_DIR}/external/cpp/rte
    ${EAM_RRTMGP_DIR}/external/cpp/rte/kernels
    ${EAM_RRTMGP_DIR}/external/cpp/rrtmgp
    ${EAM_RRTMGP_DIR}/external/cpp/rrtmgp/kernels
    ${NetCDF_C_PATH}/include
)

# Build RRTMGP interface; note that we separate the SCREAM-specific RRTMGP interface
# from the external core RRTMGP library because, ideally, the RRTMGP library has its
# own build, and we would just use add_subdirectory() above to build it. Also, this
# separates out the code that comprises the core RRTMGP library from the extensions
# and examples that we have modified for use in SCREAM specifically.

set(SCREAM_RRTMGP_SOURCES_INTERFACE
  eamxx_rrtmgp_interface.cpp
)

add_library(eamxx_rrtmgp_interface ${SCREAM_RRTMGP_SOURCES_INTERFACE})

find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATH} PATH_SUFFIXES lib lib64)
target_link_libraries(eamxx_rrtmgp_interface PUBLIC ${NETCDF_C} rrtmgp scream_share Kokkos::kokkos)
target_include_directories(eamxx_rrtmgp_interface PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(eamxx_rrtmgp_interface SYSTEM PUBLIC
    ${NetCDF_C_PATH}/include
    ${EAM_RRTMGP_DIR}/external)

##################################
#        SCREAM_RRTMGP           #
##################################

set(SCREAM_RRTMGP_SOURCES
  eamxx_rrtmgp_process_interface.cpp
  shr_orb_mod_c2f.F90
)

add_library(scream_rrtmgp ${SCREAM_RRTMGP_SOURCES})
target_link_libraries(scream_rrtmgp PUBLIC scream_share physics_share csm_share eamxx_rrtmgp_interface Kokkos::kokkos)
set_target_properties(scream_rrtmgp PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(scream_rrtmgp PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Ensure RRTMGP lookup tables are present in the data dir
set (RRTMGP_TABLES
  scream/init/rrtmgp-data-sw-g112-210809.nc
  scream/init/rrtmgp-data-lw-g128-210809.nc
  scream/init/rrtmgp-allsky.nc
  scream/init/rrtmgp-cloud-optics-coeffs-sw.nc
  scream/init/rrtmgp-cloud-optics-coeffs-lw.nc
  scream/init/rrtmgp-data-sw-g224-2018-12-04.nc
  scream/init/rrtmgp-data-lw-g256-2018-12-04.nc
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

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE scream_rrtmgp)
endif()
