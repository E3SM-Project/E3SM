module elm_instMod

  use atm2lndType       , only : atm2lnd_type
  use CanopyStateType   , only : canopystate_type
  use CNCarbonStateType , only : carbonstate_type
  use ChemStateType     , only : chemstate_type
  use decompMod         , only : bounds_type
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use ColumnDataType    , only : col_ws, col_wf, col_ef, col_es
  use VegetationDataType, only : veg_wf
  use TopounitDataType  , only : top_as
  use ColumnType        , only : col_pp
  
  implicit none
  save
  public

  type(atm2lnd_type)       :: atm2lnd_vars
  type(canopystate_type)   :: canopystate_vars
  type(carbonstate_type)   :: carbonstate_vars
  type(chemstate_type)     :: chemstate_vars
  type(energyflux_type)    :: energyflux_vars
  type(soilhydrology_type) :: soilhydrology_vars
  type(soilstate_type)     :: soilstate_vars
  type(temperature_type)   :: temperature_vars
  type(waterflux_type)     :: waterflux_vars
  type(waterstate_type)    :: waterstate_vars

  public :: elm_inst_biogeophys

  contains

  !-----------------------------------------------------------------------

    subroutine elm_inst_biogeophys(bounds_proc)

    ! !ARGUMENTS
    implicit none

    type(bounds_type), intent(in) :: bounds_proc

    call col_pp%Init(bounds_proc%begc, bounds_proc%endc)
    call col_ws%Init(bounds_proc%begc, bounds_proc%endc)
    call col_wf%Init(bounds_proc%begc, bounds_proc%endc)
    call col_ef%Init(bounds_proc%begc, bounds_proc%endc)
    call col_es%Init(bounds_proc%begc, bounds_proc%endc)
    call veg_wf%Init(bounds_proc%begc, bounds_proc%endc)
    call top_as%Init(bounds_proc%begt, bounds_proc%endt)

    call atm2lnd_vars%Init( bounds_proc )
    call canopystate_vars%init(bounds_proc)
    call carbonstate_vars%init(bounds_proc)
    call chemstate_vars%Init(bounds_proc)
    call soilstate_vars%init(bounds_proc)
    call soilhydrology_vars%Init(bounds_proc)

  end subroutine elm_inst_biogeophys

end module elm_instMod
