list(APPEND drv_sources
  component_type_mod.F90
  map_glc2lnd_mod.F90
  map_lnd2glc_mod.F90
  map_lnd2rof_irrig_mod.F90
  seq_map_mod.F90
  seq_map_type_mod.F90
  vertical_gradient_calculator_base.F90
  vertical_gradient_calculator_2nd_order.F90
  vertical_gradient_calculator_factory.F90
  )

sourcelist_to_parent(drv_sources)