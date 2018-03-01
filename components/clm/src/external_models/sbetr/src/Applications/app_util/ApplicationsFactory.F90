module ApplicationsFactory
  !
  ! !DESCRIPTION:
  !  factory to load the specific bgc reaction modules
  !
  ! History:
  !  Created by Jinyun Tang, April 28, 2016
  !
  !
  ! !USES:
  !
  use bshr_kind_mod          , only : r8 => shr_kind_r8
  use bshr_log_mod           , only : errMsg => shr_log_errMsg
  use BGCReactionsMod        , only : bgc_reaction_type
  use PlantSoilBGCMod        , only : plant_soilbgc_type

  implicit none
  character(len=*), parameter :: mod_filename = &
       __FILE__
  private
  public :: create_betr_usr_application
  public :: AppLoadParameters

contains

  subroutine create_betr_usr_application(bgc_reaction, plant_soilbgc, method, bstatus)
  !DESCRIPTION
  !create betr applications
  !
  use betr_constants  , only : betr_errmsg_len
  use BetrStatusType  , only : betr_status_type
  implicit none
  class(bgc_reaction_type),  allocatable, intent(out) :: bgc_reaction
  class(plant_soilbgc_type), allocatable, intent(out) :: plant_soilbgc
  character(len=*),                       intent(in)  :: method
  type(betr_status_type), intent(out) :: bstatus


  call create_bgc_reaction_type(bgc_reaction, method,bstatus)

  if(bstatus%check_status())return

  call create_plant_soilbgc_type(plant_soilbgc, method,bstatus)

  end subroutine create_betr_usr_application

!-------------------------------------------------------------------------------

  subroutine create_bgc_reaction_type(bgc_reaction, method, bstatus)
    !
    ! !DESCRIPTION:
    ! create and return an object of bgc_reaction
    !
    ! !USES:
    use BGCReactionsMod , only : bgc_reaction_type
    use betr_ctrl       , only : iulog  => biulog
    use betr_constants  , only : betr_errmsg_len
    use BetrStatusType  , only : betr_status_type
    use BGCReactionsCentECACnpType, only : bgc_reaction_CENTURY_ECACNP_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_type),  allocatable, intent(inout) :: bgc_reaction
    character(len=*), intent(in)          :: method
    type(betr_status_type), intent(out)   :: bstatus

    character(len=*), parameter           :: subname = 'create_bgc_reaction_type'
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    select case(trim(method))
    case ("eca_cnp")
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_ECACNP_type())
    case default
       write(msg,*)subname //' ERROR: unknown method: ', method
       msg = trim(msg)//new_line('A')//errMsg(mod_filename, __LINE__)
       call bstatus%set_msg(msg=msg,err=-1)
    end select
  end subroutine create_bgc_reaction_type
  !-------------------------------------------------------------------------------

  subroutine create_plant_soilbgc_type(plant_soilbgc, method, bstatus)

  !DESCRIPTION
  !create and return an object of plant_soilbgc
  !
  !USES
  use PlantSoilBGCMod , only : plant_soilbgc_type
  use betr_ctrl       , only : iulog  => biulog
  use betr_constants  , only : betr_errmsg_len
  use BetrStatusType  , only : betr_status_type
  use PlantSoilBgcCnpType, only : plant_soilbgc_cnp_type

  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_type), allocatable, intent(inout) :: plant_soilbgc
  character(len=*), intent(in)          :: method
  type(betr_status_type), intent(out)   :: bstatus


  character(len=*)          , parameter   :: subname = 'create_plant_soilbgc_type'
  character(len=betr_errmsg_len) :: msg

  call bstatus%reset()

  select case(trim(method))
  case ("eca_cnp")
     allocate(plant_soilbgc, source=plant_soilbgc_cnp_type())
  case default
     write(msg, *)subname //' ERROR: unknown method: ', method
     msg = trim(msg)//new_line('A')//errMsg(mod_filename, __LINE__)
     call bstatus%set_msg(msg=msg,err=-1)
  end select

  end subroutine create_plant_soilbgc_type


  !-------------------------------------------------------------------------------
  subroutine AppLoadParameters(bgc_namelist_buffer, reaction_method, bstatus)
  !
  ! DESCRIPTION
  ! read in the parameters for specified bgc implementation
  use BiogeoConType, only : bgc_con_eca
  use betr_constants , only : betr_namelist_buffer_size_ext
  use BetrStatusType , only : betr_status_type
  implicit none
  character(len=betr_namelist_buffer_size_ext), intent(in) :: bgc_namelist_buffer
  character(len=*), intent(in) :: reaction_method
  type(betr_status_type), intent(out)   :: bstatus
  character(len=255) :: msg

  call bstatus%reset()
   select case (trim(reaction_method))
   case ("eca_cnp")
     call  bgc_con_eca%Init(bgc_namelist_buffer, bstatus)
     !do nothing
   case default
     if(trim(bgc_namelist_buffer)=='none')then
       !do nothing
     else
       msg = "no parameter file to read for the specified bgc method"//errmsg(__FILE__, __LINE__)
       call bstatus%set_msg(msg=msg,err=-1)
     endif
   end select

  end subroutine  AppLoadParameters
end module ApplicationsFactory
