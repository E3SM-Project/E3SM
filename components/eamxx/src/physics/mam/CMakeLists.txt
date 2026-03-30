if (SCREAM_CIME_BUILD)
  # PROJECT_SOURCE_DIR is SCREAM_ROOT/components
  set(EXTERNALS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../externals")
else()
  # PROJECT_SOURCE_DIR is SCREAM_ROOT/components/eamxx
  set(EXTERNALS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../../externals")
endif()

# configure and build mam4xx (C++ port of MAM4)
if (SCREAM_DOUBLE_PRECISION)
  set(MAM4XX_PRECISION "double" CACHE BOOL "Enable double precision for mam4xx")
else()
  set(MAM4XX_PRECISION "single" CACHE BOOL "Enable single precision for mam4xx")
endif()
set(MAM4XX_ENABLE_GPU ${EAMXX_ENABLE_GPU} CACHE BOOL "Enable mam4xx GPU configuration" FORCE)
set(MAM4XX_ENABLE_TESTS OFF CACHE BOOL "Disable mam4xx tests within E3SM" FORCE)
set(MAM4XX_ENABLE_SKYWALKER OFF CACHE BOOL "Disable mam4xx validation tests" FORCE)
set(MAM4XX_NUM_VERTICAL_LEVELS ${SCREAM_NUM_VERTICAL_LEV})

add_subdirectory(${EXTERNALS_SOURCE_DIR}/mam4xx ${CMAKE_BINARY_DIR}/externals/mam4xx)

# EAMxx mam4xx-based atmospheric processes
add_library(mam
  readfiles/vertical_remapper_mam4.cpp
  readfiles/vertical_remapper_exo_coldens.cpp
  eamxx_mam_generic_process_interface.cpp
  eamxx_mam_microphysics_process_interface.cpp
  eamxx_mam_microphysics_process_functions.cpp
  ${SCREAM_BASE_DIR}/src/physics/rrtmgp/shr_orb_mod_c2f.F90
  eamxx_mam_optics_process_interface.cpp
  eamxx_mam_dry_deposition_process_interface.cpp
  eamxx_mam_aci_process_interface.cpp
  eamxx_mam_wetscav_process_interface.cpp
  eamxx_mam_srf_and_online_emissions_process_interface.cpp
  eamxx_mam_constituent_fluxes_interface.cpp)
target_compile_definitions(mam PUBLIC EAMXX_HAS_MAM)
target_include_directories(mam PUBLIC
  ${mam4xx_SOURCE_DIR}/mam4xx/src
  ${mam4xx_BINARY_DIR}/mam4xx/src
)
target_link_libraries(mam PUBLIC eamxx_physics_share csm_share scream_share mam4xx)

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE mam)
endif()
