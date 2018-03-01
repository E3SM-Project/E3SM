module clm_instMod
  use ColumnType         , only : col
  use PatchType          , only : pft
  use LandunitType       , only : lun
  use decompMod          , only : bounds_type
  use TemperatureType    , only : temperature_type
  use WaterstateType     , only : Waterstate_Type
  use SoilStateType      , only : soilstate_type
  use WaterfluxType      , only : waterflux_type
  use ColumnType         , only : column_type
  use ChemStateType      , only : chemstate_type
  use SoilHydrologyType  , only : soilhydrology_type
  use atm2lndType        , only : atm2lnd_type
  use CNCarbonFluxType       , only : carbonflux_type
  use CNCarbonStateType      , only : carbonstate_type
  use CNStateType            , only : cnstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use CanopyStateType        , only : canopystate_type
  use PhosphorusFluxType     , only : phosphorusflux_type
  use PhosphorusStateType    , only : phosphorusstate_type
  use PlantMicKineticsMod    , only : PlantMicKinetics_type
  use SoilWaterRetentionCurveFactoryMod, only : create_soil_water_retention_curve
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type

  implicit none
  save
  public


  type(temperature_type)      :: temperature_vars
  type(Waterstate_Type)       :: waterstate_vars   !column water state
  type(waterflux_type)        :: waterflux_vars    ! column water flux
  type(soilstate_type)        :: soilstate_vars    !column physical state variables
  type(chemstate_type)        :: chemstate_vars    !column chemical state variables
  type(soilhydrology_type)    :: soilhydrology_vars!
  type(atm2lnd_type)          :: atm2lnd_vars
  type(carbonstate_type)      :: carbonstate_vars, c13state_vars, c14state_vars
  type(carbonflux_type)       :: carbonflux_vars, c13_cflx_vars, c14_cflx_vars
  type(nitrogenstate_type)    :: nitrogenstate_vars
  type(nitrogenflux_type)     :: nitrogenflux_vars
  type(cnstate_type)          :: cnstate_vars
  type(canopystate_type)      :: canopystate_vars
  type(phosphorusstate_type)  :: phosphorusstate_vars
  type(phosphorusflux_type)   :: phosphorusflux_vars
  type(PlantMicKinetics_type) :: PlantMicKinetics_vars
  class(soil_water_retention_curve_type), allocatable :: soil_water_retention_curve
  contains

  subroutine clm_inst(bounds)
    !
    ! !DESCRIPTION:
    ! CLM initialization - first phase
  use pftvarcon, only : pftconrd
  implicit none
  type(bounds_type), intent(in) :: bounds

  call pft%Init(bounds)
  call col%Init(bounds)
  call lun%Init(bounds)

  call pftconrd

  call temperature_vars%Init(bounds)

  call waterstate_vars%Init(bounds)

  call waterflux_vars%Init(bounds)

  call soilstate_vars%Init(bounds)

  call chemstate_vars%Init(bounds)

  call atm2lnd_vars%Init(bounds)

  call soilhydrology_vars%Init(bounds)

  call carbonflux_vars%Init(bounds)

  call nitrogenstate_vars%Init(bounds)

  call nitrogenflux_vars%Init(bounds)

  call cnstate_vars%Init(bounds)

  call canopystate_vars%Init(bounds)

  call carbonstate_vars%Init(bounds)

  call phosphorusstate_vars%Init(bounds)

  call phosphorusflux_vars%Init(bounds)

  allocate(soil_water_retention_curve, &
       source=create_soil_water_retention_curve())


  end subroutine clm_inst
end module clm_instMod
