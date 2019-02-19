module clm_instMod

  use WaterstateType     , only : Waterstate_Type
  use SoilStateType      , only : soilstate_type
  use WaterfluxType      , only : waterflux_type
  use SoilHydrologyType  , only : soilhydrology_type
  use decompMod          , only : bounds_type

  implicit none
  save
  public


  type(soilstate_type)                                :: soilstate_vars
  type(soilhydrology_type)                            :: soilhydrology_vars
  type(waterflux_type)                                :: waterflux_vars
  type(waterstate_type)                               :: waterstate_vars

  public :: clm_inst_biogeophys

  contains

  !-----------------------------------------------------------------------

    subroutine clm_inst_biogeophys(bounds_proc)

    ! !ARGUMENTS
    implicit none

    type(bounds_type), intent(in) :: bounds_proc

    call soilstate_vars%init(bounds_proc)
    call waterflux_vars%init(bounds_proc)
    !call soilhydrology_vars%Init(bounds_proc, nlfilename)
    !call waterstate_vars%init(bounds_proc,         &
    !     h2osno_col(begc:endc),                    &
    !     snow_depth_col(begc:endc),                &
    !     soilstate_vars%watsat_col(begc:endc, 1:), &
    !     temperature_vars%t_soisno_col(begc:endc, -nlevsno+1:) )

  end subroutine clm_inst_biogeophys

end module clm_instMod
