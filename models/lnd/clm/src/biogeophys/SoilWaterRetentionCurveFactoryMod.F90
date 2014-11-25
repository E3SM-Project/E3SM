module SoilWaterRetentionCurveFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of soil_water_retention_curve_type. This module figures
  ! out the particular type to return.
  !
  ! !USES:
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_varctl          , only : iulog
  implicit none
  save
  private
  !
  ! !PUBLIC ROUTINES:
  public :: create_soil_water_retention_curve  ! create an object of class soil_water_retention_curve_type

contains

  !-----------------------------------------------------------------------
  function create_soil_water_retention_curve() result(soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Create and return an object of soil_water_retention_curve_type. The particular type
    ! is determined based on a namelist parameter.
    !
    ! !USES:
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SoilWaterRetentionCurveClappHornberg1978Mod, only : soil_water_retention_curve_clapp_hornberg_1978_type
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_type), allocatable :: soil_water_retention_curve  ! function result
    !
    ! !LOCAL VARIABLES:

    ! For now, hard-code the method. Eventually this will be set from namelist, either by
    ! this routine (appropriate if the 'method' is in its own namelist group), or do the
    ! namelist read outside this module and pass the method in as a parameter (appropriate
    ! if the 'method' is part of a larger namelist group).
    character(len=*), parameter :: method = "clapphornberg_1978"
    
    character(len=*), parameter :: subname = 'create_soil_water_retention_curve'
    !-----------------------------------------------------------------------
    
    select case (method)
       
    case ("clapphornberg_1978")
       allocate(soil_water_retention_curve, &
            source=soil_water_retention_curve_clapp_hornberg_1978_type())

    case default
       write(iulog,*) subname//' ERROR: unknown method: ', method
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end function create_soil_water_retention_curve

end module SoilWaterRetentionCurveFactoryMod
