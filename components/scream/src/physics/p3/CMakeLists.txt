set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(P3_SRCS
  p3_f90.cpp
  p3_functions_f90.cpp
  p3_ic_cases.cpp
  p3_iso_c.f90
  p3_iso_f.f90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/micro_p3.F90
  atmosphere_microphysics.cpp
  p3_inputs_initializer.cpp
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
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
    p3_prevent_ice_overdepletion.cpp
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
    p3_ice_supersat_conservation.cpp
    p3_nc_conservation.cpp
    p3_nr_conservation.cpp
    p3_ni_conservation.cpp
    p3_water_vapor_conservation.cpp
  ) # P3 ETI SRCS
endif()

add_library(p3 ${P3_SRCS})
target_include_directories(p3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_include_directories(p3 SYSTEM PUBLIC ${CIMEROOT}/src/share/include)
set_target_properties(p3 PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(p3 physics_share scream_share)
target_compile_options(p3 PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
