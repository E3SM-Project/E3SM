if (SCREAM_CIME_BUILD)
  # PROJECT_SOURCE_DIR is SCREAM_ROOT/components
  set(EXTERNALS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../externals")
else()
  # PROJECT_SOURCE_DIR is SCREAM_ROOT/components/eamxx
  set(EXTERNALS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../externals")
endif()

# Haero (Kokkos-based aerosol interface library)
set(HAERO_ENABLE_GPU ${EAMXX_ENABLE_GPU})
set(HAERO_ENABLE_MPI ON)
set(HAERO_ENABLE_TESTS OFF)
if (SCREAM_DOUBLE_PRECISION)
  set(HAERO_PRECISION "double")
else()
  set(HAERO_PRECISION "single")
endif()
list(APPEND CMAKE_MODULE_PATH
  ${EXTERNALS_SOURCE_DIR}/haero/cmake
)
add_subdirectory(${EXTERNALS_SOURCE_DIR}/haero ${CMAKE_BINARY_DIR}/externals/haero)

# mam4xx (C++ port of MAM4)
set(ENABLE_TESTS OFF)
set(ENABLE_SKYWALKER OFF)
set(NUM_VERTICAL_LEVELS ${SCREAM_NUM_VERTICAL_LEV})
set(HAERO_C_STANDARD ${CMAKE_C_STANDARD})
set(HAERO_C_COMPILER ${CMAKE_C_COMPILER})
set(HAERO_CXX_STANDARD ${CMAKE_CXX_STANDARD})
set(HAERO_CXX_COMPILER ${CMAKE_CXX_COMPILER})
list(APPEND CMAKE_MODULE_PATH
  ${EXTERNALS_SOURCE_DIR}/mam4xx/cmake
)
#include_directories(
  #  ${EXTERNALS_SOURCE_DIR}/ekat/src
  #  ${EXTERNALS_SOURCE_DIR}/haero
  #  ${CMAKE_BINARY_DIR}/externals/ekat/src
  #  ${CMAKE_BINARY_DIR}/externals/haero
  #)
add_subdirectory(${EXTERNALS_SOURCE_DIR}/mam4xx ${CMAKE_BINARY_DIR}/externals/mam4xx)

# EAMxx mam4xx-based atmospheric processes
add_library(mam
  eamxx_mam_microphysics_process_interface.cpp
  ${SCREAM_BASE_DIR}/src/physics/rrtmgp/shr_orb_mod_c2f.F90
  eamxx_mam_optics_process_interface.cpp
  eamxx_mam_dry_deposition_process_interface.cpp
  eamxx_mam_aci_process_interface.cpp
  eamxx_mam_wetscav_process_interface.cpp
  eamxx_mam_srf_and_online_emissions_process_interface.cpp
  eamxx_mam_constituent_fluxes_interface.cpp)
target_compile_definitions(mam PUBLIC EAMXX_HAS_MAM)
add_dependencies(mam mam4xx)
target_include_directories(mam PUBLIC
  ${EXTERNALS_SOURCE_DIR}/haero
  ${EXTERNALS_SOURCE_DIR}/mam4xx/src
)
target_link_libraries(mam PUBLIC physics_share csm_share scream_share mam4xx haero)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE mam)
endif()
