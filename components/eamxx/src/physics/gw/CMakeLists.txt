set(GW_SRCS
  eamxx_gw_process_interface.cpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND GW_SRCS
    eti/gw_gwd_compute_tendencies_from_stress_divergence.cpp
    eti/gw_gw_prof.cpp
    eti/gw_momentum_energy_conservation.cpp
    eti/gw_gwd_compute_stress_profiles_and_diffusivities.cpp
    eti/gw_gwd_project_tau.cpp
    eti/gw_gwd_precalc_rhoi.cpp
    eti/gw_gw_drag_prof.cpp
    eti/gw_gw_front_project_winds.cpp
    eti/gw_gw_front_gw_sources.cpp
    eti/gw_gw_cm_src.cpp
    eti/gw_gw_convect_project_winds.cpp
    eti/gw_gw_heating_depth.cpp
    eti/gw_gw_storm_speed.cpp
    eti/gw_gw_convect_gw_sources.cpp
    eti/gw_gw_beres_src.cpp
    eti/gw_gw_ediff.cpp
    eti/gw_gw_diff_tend.cpp
    eti/gw_gw_oro_src.cpp
    eti/gw_gw_common_init.cpp
    eti/gw_vd_lu_decomp.cpp
    eti/gw_vd_lu_solve.cpp
  ) # GW ETI SRCS
endif()

add_library(gw ${GW_SRCS})

target_include_directories(gw PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
)
target_link_libraries(gw physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE gw)
endif()
