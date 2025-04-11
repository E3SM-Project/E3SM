
# build_options.mk stuff handled here
list(APPEND CPPDEFS "-DCORE_SEAICE" "-Dcoupled" "-DCCSMCOUPLED")
list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_seaice/icepack/columnphysics" "${CMAKE_BINARY_DIR}/core_seaice/column" "${CMAKE_BINARY_DIR}/core_seaice/shared" "${CMAKE_BINARY_DIR}/core_seaice/analysis_members" "${CMAKE_BINARY_DIR}/core_seaice/model_forward")


# driver (files live in E3SM)
list(APPEND RAW_SOURCES
  ../../mpas-seaice/driver/ice_comp_mct.F
  ../../mpas-seaice/driver/mpassi_cpl_indices.F
  ../../mpas-seaice/driver/mpassi_mct_vars.F
)

# icepack
list(APPEND RAW_SOURCES
  core_seaice/icepack/columnphysics/icepack_aerosol.F90
  core_seaice/icepack/columnphysics/icepack_age.F90
  core_seaice/icepack/columnphysics/icepack_algae.F90
  core_seaice/icepack/columnphysics/icepack_atmo.F90
  core_seaice/icepack/columnphysics/icepack_brine.F90
  core_seaice/icepack/columnphysics/icepack_firstyear.F90
  core_seaice/icepack/columnphysics/icepack_flux.F90
  core_seaice/icepack/columnphysics/icepack_fsd.F90
  core_seaice/icepack/columnphysics/icepack_intfc.F90
  core_seaice/icepack/columnphysics/icepack_isotope.F90
  core_seaice/icepack/columnphysics/icepack_itd.F90
  core_seaice/icepack/columnphysics/icepack_kinds.F90
  core_seaice/icepack/columnphysics/icepack_mechred.F90
  core_seaice/icepack/columnphysics/icepack_meltpond_lvl.F90
  core_seaice/icepack/columnphysics/icepack_meltpond_topo.F90
  core_seaice/icepack/columnphysics/icepack_mushy_physics.F90
  core_seaice/icepack/columnphysics/icepack_ocean.F90
  core_seaice/icepack/columnphysics/icepack_orbital.F90
  core_seaice/icepack/columnphysics/icepack_parameters.F90
  core_seaice/icepack/columnphysics/icepack_shortwave.F90
  core_seaice/icepack/columnphysics/icepack_shortwave_data.F90
  core_seaice/icepack/columnphysics/icepack_snow.F90
  core_seaice/icepack/columnphysics/icepack_therm_bl99.F90
  core_seaice/icepack/columnphysics/icepack_therm_itd.F90
  core_seaice/icepack/columnphysics/icepack_therm_mushy.F90
  core_seaice/icepack/columnphysics/icepack_therm_shared.F90
  core_seaice/icepack/columnphysics/icepack_therm_vertical.F90
  core_seaice/icepack/columnphysics/icepack_tracers.F90
  core_seaice/icepack/columnphysics/icepack_warnings.F90
  core_seaice/icepack/columnphysics/icepack_wavefracspec.F90
  core_seaice/icepack/columnphysics/icepack_zbgc.F90
  core_seaice/icepack/columnphysics/icepack_zbgc_shared.F90
)

# column
list(APPEND RAW_SOURCES
  core_seaice/column/ice_colpkg.F90
  core_seaice/column/ice_kinds_mod.F90
  core_seaice/column/ice_warnings.F90
  core_seaice/column/ice_colpkg_shared.F90
  core_seaice/column/constants/cesm/ice_constants_colpkg.F90
  core_seaice/column/ice_therm_shared.F90
  core_seaice/column/ice_orbital.F90
  core_seaice/column/ice_mushy_physics.F90
  core_seaice/column/ice_therm_mushy.F90
  core_seaice/column/ice_atmo.F90
  core_seaice/column/ice_age.F90
  core_seaice/column/ice_firstyear.F90
  core_seaice/column/ice_flux_colpkg.F90
  core_seaice/column/ice_meltpond_cesm.F90
  core_seaice/column/ice_meltpond_lvl.F90
  core_seaice/column/ice_meltpond_topo.F90
  core_seaice/column/ice_therm_vertical.F90
  core_seaice/column/ice_therm_bl99.F90
  core_seaice/column/ice_therm_0layer.F90
  core_seaice/column/ice_itd.F90
  core_seaice/column/ice_colpkg_tracers.F90
  core_seaice/column/ice_therm_itd.F90
  core_seaice/column/ice_shortwave.F90
  core_seaice/column/ice_mechred.F90
  core_seaice/column/ice_aerosol.F90
  core_seaice/column/ice_brine.F90
  core_seaice/column/ice_algae.F90
  core_seaice/column/ice_zbgc.F90
  core_seaice/column/ice_zbgc_shared.F90
  core_seaice/column/ice_zsalinity.F90
  core_seaice/column/ice_snow.F90
)

# shared
list(APPEND RAW_SOURCES
  core_seaice/shared/mpas_seaice_time_integration.F
  core_seaice/shared/mpas_seaice_advection_incremental_remap_tracers.F
  core_seaice/shared/mpas_seaice_advection_incremental_remap.F
  core_seaice/shared/mpas_seaice_advection_upwind.F
  core_seaice/shared/mpas_seaice_advection.F
  core_seaice/shared/mpas_seaice_velocity_solver.F
  core_seaice/shared/mpas_seaice_velocity_solver_weak.F
  core_seaice/shared/mpas_seaice_velocity_solver_variational.F
  core_seaice/shared/mpas_seaice_velocity_solver_wachspress.F
  core_seaice/shared/mpas_seaice_velocity_solver_pwl.F
  core_seaice/shared/mpas_seaice_velocity_solver_variational_shared.F
  core_seaice/shared/mpas_seaice_velocity_solver_constitutive_relation.F
  core_seaice/shared/mpas_seaice_triangle_quadrature.F
  core_seaice/shared/mpas_seaice_wachspress_basis.F
  core_seaice/shared/mpas_seaice_forcing.F
  core_seaice/shared/mpas_seaice_initialize.F
  core_seaice/shared/mpas_seaice_testing.F
  core_seaice/shared/mpas_seaice_mesh.F
  core_seaice/shared/mpas_seaice_diagnostics.F
  core_seaice/shared/mpas_seaice_numerics.F
  core_seaice/shared/mpas_seaice_constants.F
  core_seaice/shared/mpas_seaice_column.F
  core_seaice/shared/mpas_seaice_icepack.F
  core_seaice/shared/mpas_seaice_diagnostics.F
  core_seaice/shared/mpas_seaice_error.F
  core_seaice/shared/mpas_seaice_mesh_pool.F
  core_seaice/shared/mpas_seaice_prescribed.F
  core_seaice/shared/mpas_seaice_special_boundaries.F
)

# analysis members
list(APPEND RAW_SOURCES
  core_seaice/analysis_members/mpas_seaice_analysis_driver.F
  core_seaice/analysis_members/mpas_seaice_high_frequency_output.F
  core_seaice/analysis_members/mpas_seaice_temperatures.F
  core_seaice/analysis_members/mpas_seaice_thicknesses.F
  core_seaice/analysis_members/mpas_seaice_regional_statistics.F
  core_seaice/analysis_members/mpas_seaice_ridging_diagnostics.F
  core_seaice/analysis_members/mpas_seaice_conservation_check.F
  core_seaice/analysis_members/mpas_seaice_geographical_vectors.F
  core_seaice/analysis_members/mpas_seaice_ice_present.F
  core_seaice/analysis_members/mpas_seaice_time_series_stats.F
  core_seaice/analysis_members/mpas_seaice_load_balance.F
  core_seaice/analysis_members/mpas_seaice_maximum_ice_presence.F
  core_seaice/analysis_members/mpas_seaice_miscellaneous.F
  core_seaice/analysis_members/mpas_seaice_area_variables.F
  core_seaice/analysis_members/mpas_seaice_pond_diagnostics.F
  core_seaice/analysis_members/mpas_seaice_deactivate_unneeded_fields.F
  core_seaice/analysis_members/mpas_seaice_pointwise_stats.F
  core_seaice/analysis_members/mpas_seaice_unit_conversion.F
  core_seaice/analysis_members/mpas_seaice_ice_shelves.F
)

# model_forward (DISABLE qsmp for these)
set(SEAICE_MODEL_FORWARD
  core_seaice/model_forward/mpas_seaice_core.F
  core_seaice/model_forward/mpas_seaice_core_interface.F
)
list(APPEND RAW_SOURCES ${SEAICE_MODEL_FORWARD})
list(APPEND DISABLE_QSMP ${SEAICE_MODEL_FORWARD})
list(APPEND NOOPT_FILES "core_seaice/icepack/columnphysics/icepack_shortwave_data.F90")

# Generate core input
handle_st_nl_gen("namelist.seaice" "streams.seaice stream_list.seaice. listed" ${CORE_INPUT_DIR} ${CORE_BLDDIR})
