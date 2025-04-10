set(PATH_TO_LEGACY_GW ${SCREAM_BASE_DIR}/../eam/src/physics/cam/gw)
set(GW_SRCS
  ${PATH_TO_LEGACY_GW}/gw_utils.F90
  ${PATH_TO_LEGACY_GW}/gw_common.F90
  ${PATH_TO_LEGACY_GW}/gw_convect.F90
  ${PATH_TO_LEGACY_GW}/gw_diffusion.F90
  ${PATH_TO_LEGACY_GW}/gw_front.F90
  ${PATH_TO_LEGACY_GW}/gw_oro.F90
  gw_iso_c.f90
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND GW_SRCS
  ) # GW ETI SRCS
endif()

add_library(gw ${GW_SRCS})
set_target_properties(gw PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(gw PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
  ${PATH_TO_LEGACY_GW}
)
target_link_libraries(gw physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE gw)
endif()
