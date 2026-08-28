set(ZM_SRCS
  # ----------------------------------------------------------------------------
  # C++/EAMxx atomsphere process for ZM
  eamxx_zm_process_interface.cpp
  # ----------------------------------------------------------------------------
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND ZM_SRCS
    eti/zm_opts.cpp
    eti/zm_input_state.cpp
    eti/zm_output_tend.cpp
    eti/zm_invert_entropy.cpp
    eti/zm_entropy.cpp
    eti/zm_transport_tracer.cpp
    eti/zm_transport_momentum.cpp
    eti/zm_compute_dilute_cape.cpp
    eti/zm_find_mse_max.cpp
    eti/zm_compute_dilute_parcel.cpp
    eti/zm_compute_cape_from_parcel.cpp
    eti/zm_conv_mcsp_calculate_shear.cpp
    eti/zm_conv_mcsp_tend.cpp
    eti/zm_conv_main.cpp
    eti/zm_conv_evap.cpp
    eti/zm_calc_fractional_entrainment.cpp
    eti/zm_downdraft_properties.cpp
    eti/zm_cloud_properties.cpp
    eti/zm_closure.cpp
    eti/zm_calc_output_tend.cpp
  ) # ZM ETI SRCS
endif()

add_library(zm ${ZM_SRCS})

set_target_properties(zm PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_compile_definitions(zm PUBLIC EAMXX_HAS_ZM)

target_include_directories(zm PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${PATH_TO_LEGACY_ZM}
  ${PATH_TO_LEGACY_ZM}/..
)
target_link_libraries(zm eamxx_physics_share scream_share csm_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE zm)
endif()
