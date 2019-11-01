
# driver (files live in E3SM)
list(APPEND RAW_SOURCES
  ../../mpas-albany-landice/driver/glc_comp_mct.F
  ../../mpas-albany-landice/driver/glc_cpl_indices.F
  ../../mpas-albany-landice/driver/glc_mct_vars.F
)

# shared
list(APPEND RAW_SOURCES
  core_landice/shared/mpas_li_constants.F
  core_landice/shared/mpas_li_mask.F
  core_landice/shared/mpas_li_setup.F
)

# analysis members
list(APPEND RAW_SOURCES
  core_landice/analysis_members/mpas_li_analysis_driver.F
  core_landice/analysis_members/mpas_li_global_stats.F
  core_landice/analysis_members/mpas_li_regional_stats.F
)

# mode forward
list(APPEND RAW_SOURCES
  core_landice/mode_forward/mpas_li_core.F
  core_landice/mode_forward/mpas_li_core_interface.F
  core_landice/mode_forward/mpas_li_time_integration.F
  core_landice/mode_forward/mpas_li_time_integration_fe.F
  core_landice/mode_forward/mpas_li_diagnostic_vars.F
  core_landice/mode_forward/mpas_li_advection.F
  core_landice/mode_forward/mpas_li_calving.F
  core_landice/mode_forward/mpas_li_statistics.F
  core_landice/mode_forward/mpas_li_velocity.F
  core_landice/mode_forward/mpas_li_thermal.F
  core_landice/mode_forward/mpas_li_iceshelf_melt.F
  core_landice/mode_forward/mpas_li_sia.F
  core_landice/mode_forward/mpas_li_velocity_simple.F
  core_landice/mode_forward/mpas_li_velocity_external.F
  core_landice/mode_forward/mpas_li_subglacial_hydro.F
)

if (CPPFLAGS MATCHES ".*MPAS_LI_BUILD_INTERFACE.*")
  list(APPEND RAW_SOURCES core_landice/mode_forward/Interface_velocity_solver.cpp)
endif()

# Generate core input
handle_st_nl_gen("namelist.landice" "streams.landice stream_list.landice. listed" ${CORE_INPUT_DIR} ${CORE_BLDDIR})
