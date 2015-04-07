# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources 
  DaylengthMod.F90
  HumanIndexMod.F90
  IrrigationMod.F90
  RootBiophysMod.F90
  SoilStateType.F90
  SoilWaterRetentionCurveMod.F90
  TemperatureType.F90
  WaterfluxType.F90
  )

sourcelist_to_parent(clm_sources)
