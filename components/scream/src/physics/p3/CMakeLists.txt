set(P3_SRCS
  p3_f90.cpp
  p3_functions_f90.cpp
  p3_ic_cases.cpp
  p3_iso_c.f90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/micro_p3.F90
  atmosphere_microphysics.cpp
  atmosphere_microphysics_run.cpp
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  list(APPEND P3_SRCS
    p3_upwind.cpp
    p3_cloud_sed.cpp
    p3_ice_sed.cpp
    p3_ice_collection.cpp
    p3_ice_melting.cpp
    p3_rain_sed.cpp
    p3_table3.cpp
    p3_table_ice.cpp
    p3_dsd2.cpp
    p3_find.cpp
    p3_update_prognostics.cpp
    p3_get_time_space_phys_variables.cpp
    p3_autoconversion.cpp
    p3_ice_deposition_sublimation.cpp
    p3_conservation.cpp
    p3_rain_self_collection.cpp
    p3_back_to_cell_average.cpp
    p3_cloud_rain_acc.cpp
    p3_calc_rime_density.cpp
    p3_cldliq_imm_freezing.cpp
    p3_rain_imm_freezing.cpp
    p3_droplet_self_coll.cpp
    p3_update_prognostics.cpp
    p3_evaporate_rain.cpp
    p3_impose_max_total_ni.cpp
    p3_calc_liq_relaxation_timescale.cpp
    p3_ice_relaxation_timescale.cpp
    p3_ice_nucleation.cpp
    p3_ice_cldliq_wet_growth.cpp
    p3_get_latent_heat.cpp
    p3_check_values.cpp
    p3_incloud_mixingratios.cpp
    p3_subgrid_variance_scaling.cpp
    p3_main.cpp
    p3_main_part1.cpp
    p3_main_part2.cpp
    p3_main_part3.cpp
    p3_ice_supersat_conservation.cpp
    p3_nc_conservation.cpp
    p3_nr_conservation.cpp
    p3_ni_conservation.cpp
    p3_prevent_liq_supersaturation.cpp
  ) # P3 ETI SRCS
endif()

add_library(p3 ${P3_SRCS})
target_include_directories(p3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
set_target_properties(p3 PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(p3 physics_share scream_share)
target_compile_options(p3 PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)

# Ensure tables are present in the data dir
if (SCREAM_DOUBLE_PRECISION)
  set (PRECISION_SUFFIX "8")
else()
  set (PRECISION_SUFFIX "4")
endif()
set (P3_TABLES
  tables/p3_lookup_table_1.dat-v4.1.1
  tables/mu_r_table_vals.dat${PRECISION_SUFFIX}
  tables/revap_table_vals.dat${PRECISION_SUFFIX}
  tables/vm_table_vals.dat${PRECISION_SUFFIX}
  tables/vn_table_vals.dat${PRECISION_SUFFIX}
)

include (ScreamUtils)
foreach (file IN ITEMS ${P3_TABLES})
  GetInputFile(${file})
endforeach()

# This executable can be used to re-generate tables in ${SCREAM_DATA_DIR}
add_executable(p3_tables_setup EXCLUDE_FROM_ALL p3_tables_setup.cpp)
target_link_libraries(p3_tables_setup p3)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
