set(REDOPT
  ../driver-mct/main/seq_io_mod.F90
  elm/src/biogeophys/BandDiagonalMod.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS REDOPT)
    e3sm_add_flags("${ITEM}" "-O1 -g")
  endforeach()
endif()

set(CICE_F90
  ice_FY.F90
  ice_aerosol.F90
  ice_age.F90
  ice_atmo.F90
  ice_blocks.F90
  ice_calendar.F90
  ice_diagnostics.F90
  ice_distribution.F90
  ice_domain.F90
  ice_domain_size.F90
  ice_dyn_evp.F90
  ice_fileunits.F90
  ice_flux.F90
  ice_forcing.F90
  ice_grid.F90
  ice_history.F90
  ice_history_fields.F90
  ice_init.F90
  ice_itd.F90
  ice_kinds_mod.F90
  ice_lvl.F90
  ice_mechred.F90
  ice_meltpond.F90
  ice_ocean.F90
  ice_orbital.F90
  ice_probability.F90
  ice_probability_tools.F90
  ice_read_write.F90
  ice_restoring.F90
  ice_shortwave.F90
  ice_spacecurve.F90
  ice_state.F90
  ice_step_mod.F90
  ice_therm_itd.F90
  ice_therm_vertical.F90
  ice_transport_driver.F90
  ice_transport_remap.F90
  ice_work.F90)

foreach(ITEM IN LISTS CICE_F90)
  e3sm_add_flags("cice/src/source/${ITEM}" "-O0")
endforeach()
