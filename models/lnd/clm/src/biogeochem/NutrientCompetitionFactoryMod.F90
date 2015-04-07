module NutrientCompetitionFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of nutrient_competition_method_type. This module figures
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
  public :: create_nutrient_competition_method  ! create an object of class nutrient_competition_method_type

contains

  !-----------------------------------------------------------------------
  function create_nutrient_competition_method() result(nutrient_competition_method)
    !
    ! !DESCRIPTION:
    ! Create and return an object of nutrient_competition_method_type. The particular type
    ! is determined based on a namelist parameter.
    !
    ! !USES:
    use NutrientCompetitionMethodMod, only : nutrient_competition_method_type
    use NutrientCompetitionCLM45defaultMod, only : nutrient_competition_clm45default_type
    !
    ! !ARGUMENTS:
    class(nutrient_competition_method_type), allocatable :: nutrient_competition_method  ! function result
    !
    ! !LOCAL VARIABLES:

    ! For now, hard-code the method. Eventually this will be set from namelist, either by
    ! this routine (appropriate if the 'method' is in its own namelist group), or do the
    ! namelist read outside this module and pass the method in as a parameter (appropriate
    ! if the 'method' is part of a larger namelist group).
    character(len=*), parameter :: method = "clm45default"
    
    character(len=*), parameter :: subname = 'create_nutrient_competition_method'
    !-----------------------------------------------------------------------
    
    select case (method)
       
    case ("clm45default")
       allocate(nutrient_competition_method, &
            source=nutrient_competition_clm45default_type())

    case default
       write(iulog,*) subname//' ERROR: unknown method: ', method
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end function create_nutrient_competition_method

end module NutrientCompetitionFactoryMod
