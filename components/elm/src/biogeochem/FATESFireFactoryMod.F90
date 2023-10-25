module FATESFireFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of fire_method_type. This module figures
  ! out the particular type to return.
  !
  ! !USES:
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use elm_varctl          , only : iulog

  implicit none
  save
  private
  !
  ! !PUBLIC ROUTINES:
  public :: create_fates_fire_data_method  ! create an object of class fates_fire_base_type

  ! These parameters set the ranges of the cases in subroutine
  ! create_fates_fire_data_method. We declare them public in order to
  ! use them as flags elsewhere in the CTSM and FATES-SPITFIRE.
  ! They correspond one-to-one to the fates_spitfire_mode options listed
  ! in namelist definition file
  integer, public, parameter :: no_fire = 0              ! value of no_fire mode
  integer, public, parameter :: scalar_lightning = 1     ! value of scalar_lightning mode
  integer, public, parameter :: lightning_from_data = 2  ! value of lightning_from_data mode
  integer, public, parameter :: successful_ignitions = 3 ! value of successful_ignitions mode
  integer, public, parameter :: anthro_ignitions = 4     ! value of anthro_ignitions mode
  integer, public, parameter :: anthro_suppression = 5   ! value of anthro_suppression mode

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine create_fates_fire_data_method( fates_fire_data_method )
    !
    ! !DESCRIPTION:
    ! Create and return an object of fates_fire_data_method_type.
    ! The particular type is determined based on a namelist parameter.
    !
    ! !USES:
    use elm_varctl,           only: fates_spitfire_mode
    use FATESFireBase,        only: fates_fire_base_type
    use FATESFireNoDataMod,   only: fates_fire_no_data_type
    use FATESFireDataMod,     only: fates_fire_data_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type), allocatable, intent(inout) :: fates_fire_data_method  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: current_case

    character(len=*), parameter :: subname = 'create_fates_fire_data_method'
    !-----------------------------------------------------------------------

    current_case = fates_spitfire_mode

    select case (current_case)

    case (no_fire:scalar_lightning)
       allocate(fates_fire_no_data_type :: fates_fire_data_method)
    case (lightning_from_data:anthro_suppression)
       allocate(fates_fire_data_type :: fates_fire_data_method)

    case default
       write(iulog,*) subname//' ERROR: unknown method: ', fates_spitfire_mode
       call endrun(msg=errMsg(sourcefile, __LINE__))

    end select

  end subroutine create_fates_fire_data_method

end module FATESFireFactoryMod
