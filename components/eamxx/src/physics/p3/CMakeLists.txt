set(P3_SRCS
  p3_f90.cpp
  p3_ic_cases.cpp
  p3_iso_c.f90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/micro_p3.F90
  atmosphere_microphysics.cpp
  atmosphere_microphysics_run.cpp
)

if (NOT SCREAM_LIB_ONLY)
  list(APPEND P3_SRCS
    p3_functions_f90.cpp
    p3_main_wrap.cpp
    ) # Add f90 bridges needed for testing
endif()

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND P3_SRCS
    eti/p3_upwind.cpp
    eti/p3_cloud_sed.cpp
    eti/p3_ice_sed.cpp
    eti/p3_ice_collection.cpp
    eti/p3_ice_melting.cpp
    eti/p3_rain_sed.cpp
    eti/p3_table3.cpp
    eti/p3_table_ice.cpp
    eti/p3_dsd2.cpp
    eti/p3_find.cpp
    eti/p3_update_prognostics.cpp
    eti/p3_get_time_space_phys_variables.cpp
    eti/p3_autoconversion.cpp
    eti/p3_ice_deposition_sublimation.cpp
    eti/p3_conservation.cpp
    eti/p3_rain_self_collection.cpp
    eti/p3_back_to_cell_average.cpp
    eti/p3_cloud_rain_acc.cpp
    eti/p3_calc_rime_density.cpp
    eti/p3_cldliq_imm_freezing.cpp
    eti/p3_rain_imm_freezing.cpp
    eti/p3_droplet_self_coll.cpp
    eti/p3_evaporate_rain.cpp
    eti/p3_impose_max_total_ni.cpp
    eti/p3_calc_liq_relaxation_timescale.cpp
    eti/p3_ice_relaxation_timescale.cpp
    eti/p3_ice_nucleation.cpp
    eti/p3_ice_cldliq_wet_growth.cpp
    eti/p3_get_latent_heat.cpp
    eti/p3_check_values.cpp
    eti/p3_incloud_mixingratios.cpp
    eti/p3_subgrid_variance_scaling.cpp
    eti/p3_main.cpp
    eti/p3_main_part1.cpp
    eti/p3_main_part2.cpp
    eti/p3_main_part3.cpp
    eti/p3_ice_supersat_conservation.cpp
    eti/p3_nc_conservation.cpp
    eti/p3_nr_conservation.cpp
    eti/p3_ni_conservation.cpp
    eti/p3_prevent_liq_supersaturation.cpp
  ) # P3 ETI SRCS
endif()

add_library(p3 ${P3_SRCS})
set_target_properties(p3 PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(p3 PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../share
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
)
target_link_libraries(p3 physics_share scream_share)

# Ensure tables are present in the data dir
if (SCREAM_DOUBLE_PRECISION)
  set (PRECISION_SUFFIX "8")
else()
  set (PRECISION_SUFFIX "4")
endif()
set (P3_TABLES
  scream/tables/p3_lookup_table_1.dat-v4.1.1
  scream/tables/mu_r_table_vals.dat${PRECISION_SUFFIX}
  scream/tables/revap_table_vals.dat${PRECISION_SUFFIX}
  scream/tables/vm_table_vals.dat${PRECISION_SUFFIX}
  scream/tables/vn_table_vals.dat${PRECISION_SUFFIX}
)

include (ScreamUtils)
foreach (file IN ITEMS ${P3_TABLES})
  GetInputFile(${file})
endforeach()

# This executable can be used to re-generate tables in ${SCREAM_DATA_DIR}
add_executable(p3_tables_setup EXCLUDE_FROM_ALL p3_tables_setup.cpp)
target_link_libraries(p3_tables_setup p3)

#crusher change
if (Kokkos_ENABLE_HIP)
set_source_files_properties(p3_functions_f90.cpp  PROPERTIES COMPILE_FLAGS -O0)
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
