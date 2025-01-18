
# build_options.mk stuff handled here
list(APPEND CPPDEFS "-DCORE_LANDICE")
list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_landice/shared" "${CMAKE_BINARY_DIR}/core_landice/analysis_members" "${CMAKE_BINARY_DIR}/core_landice/mode_forward")

list(APPEND LIBRARIES ${Albany_LIBRARIES})

#
# Check if building with LifeV, Albany, and/or PHG external libraries
#

if (LIFEV)
  # LifeV can solve L1L2 or FO
  list(APPEND CPPDEFS "-DLIFEV" "-DUSE_EXTERNAL_L1L2" "-DUSE_EXTERNAL_FIRSTORDER" "-DMPAS_LI_BUILD_INTERFACE")
endif()

# Albany can only solve FO at present
if (ALBANY)
  list(APPEND CPPDEFS "-DUSE_EXTERNAL_FIRSTORDER" "-DMPAS_LI_BUILD_INTERFACE")
endif()

if (LIFEV AND ALBANY)
  message(FATAL "Compiling with both LifeV and Albany is not allowed at this time.")
endif()

# PHG currently requires LifeV
if (PHG AND NOT LIFEV)
  message(FATAL "Compiling with PHG requires LifeV at this time.")
endif()

# PHG can only Stokes at present
if (PHG)
  list(APPEND CPPDEFS "-DUSE_EXTERNAL_STOKES" "-DMPAS_LI_BUILD_INTERFACE")
endif()

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
  core_landice/shared/mpas_li_mesh.F
  core_landice/shared/mpas_li_config.F
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
  core_landice/mode_forward/mpas_li_time_integration_fe_rk.F
  core_landice/mode_forward/mpas_li_diagnostic_vars.F
  core_landice/mode_forward/mpas_li_advection.F
  core_landice/mode_forward/mpas_li_advection_fct_shared.F
  core_landice/mode_forward/mpas_li_advection_fct.F
  core_landice/mode_forward/mpas_li_calving.F
  core_landice/mode_forward/mpas_li_statistics.F
  core_landice/mode_forward/mpas_li_velocity.F
  core_landice/mode_forward/mpas_li_thermal.F
  core_landice/mode_forward/mpas_li_iceshelf_melt.F
  core_landice/mode_forward/mpas_li_sia.F
  core_landice/mode_forward/mpas_li_velocity_simple.F
  core_landice/mode_forward/mpas_li_velocity_external.F
  core_landice/mode_forward/mpas_li_subglacial_hydro.F
  core_landice/mode_forward/mpas_li_bedtopo.F
  core_landice/mode_forward/mpas_li_ocean_extrap.F
)

if (CPPDEFS MATCHES ".*MPAS_LI_BUILD_INTERFACE.*")
  list(APPEND RAW_SOURCES core_landice/mode_forward/Interface_velocity_solver.cpp)
endif()

# Generate core input
handle_st_nl_gen("namelist.landice" "streams.landice stream_list.landice. listed" ${CORE_INPUT_DIR} ${CORE_BLDDIR})
