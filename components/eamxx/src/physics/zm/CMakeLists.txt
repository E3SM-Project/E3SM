set(PATH_TO_LEGACY_ZM ${SCREAM_BASE_DIR}/../eam/src/physics/cam/)

set(ZM_F90_SRCS
  # ----------------------------------------------------------------------------
  # EAMxx side bridge modules
  ${CMAKE_CURRENT_SOURCE_DIR}/zm_eamxx_bridge_main.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/zm_eamxx_bridge_params.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/zm_eamxx_bridge_methods.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/zm_eamxx_bridge_physconst.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/zm_eamxx_bridge_wv_saturation.F90
  # ----------------------------------------------------------------------------
  # Fortran ZM code
  ${PATH_TO_LEGACY_ZM}/zm_conv.F90
  ${PATH_TO_LEGACY_ZM}/zm_conv_cape.F90
  ${PATH_TO_LEGACY_ZM}/zm_conv_util.F90
  ${PATH_TO_LEGACY_ZM}/zm_conv_types.F90
  ${PATH_TO_LEGACY_ZM}/zm_aero_type.F90
  ${PATH_TO_LEGACY_ZM}/zm_microphysics_state.F90
  # ----------------------------------------------------------------------------
  # shared fotran modules
  ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in
  # ----------------------------------------------------------------------------
  # misc other fotran modules
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/wv_sat_methods.F90
  ${SCREAM_BASE_DIR}/../eam/src/utils/cam_abortutils.F90
  ${SCREAM_BASE_DIR}/../eam/src/control/cam_logfile.F90
  # ----------------------------------------------------------------------------
)

set(ZM_CXX_SRCS
  # ----------------------------------------------------------------------------
  # C++/EAMxx atomsphere process for ZM
  eamxx_zm_process_interface.cpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND ZM_SRCS
  )
endif()

add_library(zm ${ZM_F90_SRCS} ${ZM_CXX_SRCS})

set_target_properties(zm PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(zm PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
  ${PATH_TO_LEGACY_ZM}
)
target_link_libraries(zm physics_share scream_share csm_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE zm)
endif()
