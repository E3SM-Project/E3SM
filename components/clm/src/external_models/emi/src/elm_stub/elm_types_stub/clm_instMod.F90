module clm_instMod

  use atm2lndType       , only : atm2lnd_type
  use CanopyStateType   , only : canopystate_type
  use ChemStateType     , only : chemstate_type
  use decompMod         , only : bounds_type
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type

  implicit none
  save
  public

  type(atm2lnd_type)       :: atm2lnd_vars
  type(canopystate_type)   :: canopystate_vars
  type(chemstate_type)     :: chemstate_vars
  type(energyflux_type)    :: energyflux_vars
  type(soilhydrology_type) :: soilhydrology_vars
  type(soilstate_type)     :: soilstate_vars
  type(temperature_type)   :: temperature_vars
  type(waterflux_type)     :: waterflux_vars
  type(waterstate_type)    :: waterstate_vars

  public :: clm_inst_biogeophys

  contains

  !-----------------------------------------------------------------------

    subroutine clm_inst_biogeophys(bounds_proc)

    ! !ARGUMENTS
    implicit none

    type(bounds_type), intent(in) :: bounds_proc

    call atm2lnd_vars%Init( bounds_proc )
    call canopystate_vars%init(bounds_proc)
    call chemstate_vars%Init(bounds_proc)
    call soilstate_vars%init(bounds_proc)
    call soilhydrology_vars%Init(bounds_proc)
    call temperature_vars%init(bounds_proc)
    call energyflux_vars%init(bounds_proc)
    call waterflux_vars%init(bounds_proc)
    call waterstate_vars%init(bounds_proc)

  end subroutine clm_inst_biogeophys

end module clm_instMod
