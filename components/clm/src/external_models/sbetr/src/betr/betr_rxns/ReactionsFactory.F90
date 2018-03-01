module ReactionsFactory
  !
  ! !DESCRIPTION:
  !  factory to load the specific bgc reaction modules
  !
  ! History:
  !  Created by Jinyun Tang, Oct 2, 2014
  !
  !
  ! !USES:
  !
  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  use BGCReactionsMod , only : bgc_reaction_type
  use PlantSoilBGCMod , only : plant_soilbgc_type
  implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  private

  public :: create_betr_def_application

contains

  subroutine create_betr_def_application(bgc_reaction, plant_soilbgc, method, yesno)
  !DESCRIPTION
  !create betr applications
  !
  implicit none
  !arguments
  class(bgc_reaction_type),  allocatable, intent(out) :: bgc_reaction
  class(plant_soilbgc_type), allocatable, intent(out) :: plant_soilbgc
  character(len=*),                       intent(in)  :: method

  logical, intent(out) :: yesno
  yesno = is_reaction_exist(method)
  if(yesno)then
    call create_bgc_reaction_type(bgc_reaction, method)
    call create_plant_soilbgc_type(plant_soilbgc, method)
  endif

  end subroutine create_betr_def_application
!-------------------------------------------------------------------------------

  subroutine create_bgc_reaction_type(bgc_reaction, method)
    !
    ! !DESCRIPTION:
    ! create and return an object of bgc_reaction
    !
    ! !USES:
    use BGCReactionsMod            , only : bgc_reaction_type
    use MockBGCReactionsType       , only : bgc_reaction_mock_run_type
    use H2OIsotopeBGCReactionsType , only : bgc_reaction_h2oiso_type
    use betr_ctrl                  , only : iulog  => biulog
    use BetrStatusType             , only : betr_status_type
    use betr_constants             , only : betr_errmsg_len
    use DIOCBGCReactionsType       , only : bgc_reaction_dioc_run_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_type),  allocatable, intent(inout) :: bgc_reaction
    character(len=*)         , intent(in)  :: method
    !local variables
    character(len=*)         , parameter   :: subname = 'create_bgc_reaction_type'

    select case(trim(method))
    case ("mock_run")
       allocate(bgc_reaction, source=bgc_reaction_mock_run_type())
    case ("h2oiso")
       allocate(bgc_reaction, source=bgc_reaction_h2oiso_type())
    case ("doc_dic")
       allocate(bgc_reaction, source=bgc_reaction_dioc_run_type())
    case default
!x       write(iulog,*)subname //' ERROR: unknown method: ', method
!x       call endrun(msg=errMsg(mod_filename, __LINE__))
       !do nothing
    end select
  end subroutine create_bgc_reaction_type
  !-------------------------------------------------------------------------------

  subroutine create_plant_soilbgc_type(plant_soilbgc, method)
  !DESCRIPTION
  !create plant soil bgc type
  !
  !USES
  use PlantSoilBGCMod            , only : plant_soilbgc_type
  use MockPlantSoilBGCType       , only : plant_soilbgc_mock_run_type
  use H2OIsotopePlantSoilBGCType , only : plant_soilbgc_h2oiso_run_type
  use betr_ctrl                  , only : iulog  => biulog
  use DIOCPlantSoilBGCType       , only : plant_soilbgc_dioc_run_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_type), allocatable, intent(inout) :: plant_soilbgc
  character(len=*)          , intent(in)  :: method
  character(len=*)          , parameter   :: subname = 'create_standalone_plant_soilbgc_type'

  select case(trim(method))
  case ("mock_run")
     allocate(plant_soilbgc, source=plant_soilbgc_mock_run_type())
  case ("h2oiso")
     allocate(plant_soilbgc, source=plant_soilbgc_h2oiso_run_type())
  case ("doc_dic")
     allocate(plant_soilbgc, source=plant_soilbgc_dioc_run_type())

  case default
!x     write(*, *)subname //' ERROR: unknown method: ', method
!x     call endrun(msg=errMsg(mod_filename, __LINE__))
   !do nothing
  end select

  end subroutine create_plant_soilbgc_type
  !-------------------------------------------------------------------------------
  function is_reaction_exist(method)result(yesno)
  !DESCRIPTION
  !determine if it is a default betr application
  use betr_ctrl                  , only : iulog  => biulog
  implicit none
  character(len=*), intent(in) :: method
  character(len=*), parameter  :: subname = 'is_reaction_exist'
  !local variable
  logical :: yesno
  select case(trim(method))
  case ("mock_run")
     yesno = .true.
  case ("h2oiso")
     yesno = .true.
  case ("doc_dic")
     yesno = .true.
  case default
     !write(iulog, *)subname //' Warning: unknown default method: '//trim(method)//'. Looking for user defined method.'
     yesno = .false.
  end select

  end function is_reaction_exist
end module ReactionsFactory
