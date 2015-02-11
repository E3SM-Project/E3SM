module BGCReactionsFactoryMod
!
! DESCRIPTION
! factory to load the specific bgc reaction modules
!
! History:
!  Created by Jinyun Tang, Oct 2, 2014
!
  !
  ! USES
  !
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type
  use shr_log_mod           , only : errMsg => shr_log_errMsg
implicit none
  save
  private  
  public ::  ctreate_bgc_reaction_type
  public ::  is_active_betr_bgc
  
  logical, private :: active_bgc   !when it is active, betr will handel the subsurface bgc
  contains 
  
  function ctreate_bgc_reaction_type(method) result(bgc_reaction)
  !
  ! ! DESCRIPTION
  ! create and return an object of bgc_reaction
  !
  ! USES
  use BGCReactionsMod           , only : bgc_reaction_type
  use BGCReactionsMockRunType   , only : bgc_reaction_mock_run_type
  use BGCReactionsCenturyType   , only : bgc_reaction_CENTURY_type
!  use BGCReactionsO18IsotopeType, only : bgc_reaction_O18ISO_type
  use abortutils                , only : endrun
  use clm_varctl                , only : iulog
  ! !ARGUMENTS
  
  class(bgc_reaction_type), allocatable :: bgc_reaction
  
  character(len=*), intent(in) :: method
  character(len=*), parameter :: subname = 'ctreate_bgc_reaction_type'
  
  active_bgc = .false.
  
  select case(trim(method))
  case ("mock_run")
    allocate(bgc_reaction, source=bgc_reaction_mock_run_type())
  case ("century_bgc")
    active_bgc = .true.
    allocate(bgc_reaction, source=bgc_reaction_CENTURY_type())
  !case ("o18_istope")    ! on hold
  !  allocate(bgc_reaction, source=bgc_reaction_O18ISO_type())
  case default
    write(iulog,*)subname //' ERROR: unknown method: ', method
    call endrun(msg=errMsg(__FILE__, __LINE__))
  end select
  end function ctreate_bgc_reaction_type
  
  
  function is_active_betr_bgc()result(ans)
  implicit none
  
  logical :: ans
  
  ans = active_bgc
  
  end function is_active_betr_bgc  
end module BGCReactionsFactoryMod
