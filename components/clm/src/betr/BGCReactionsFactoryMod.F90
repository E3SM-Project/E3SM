module BGCReactionsFactoryMod
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
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  implicit none
  save
  private
  public ::  ctreate_bgc_reaction_type


contains

  function ctreate_bgc_reaction_type(method) result(bgc_reaction)
    !
    ! !DESCRIPTION:
    ! create and return an object of bgc_reaction
    !
    ! !USES:
    use BGCReactionsMod             , only : bgc_reaction_type
    use BGCReactionsMockRunType     , only : bgc_reaction_mock_run_type
    use BGCReactionsCenturyType     , only : bgc_reaction_CENTURY_type
    use BGCReactionsCenturyCLMType  , only : bgc_reaction_CENTURY_clm_type
    use BGCReactionsCenturyECAType  , only : bgc_reaction_CENTURY_ECA_type
    use BGCReactionsSminNType       , only : bgc_reaction_sminn_type
    use BGCReactionsCenturyCLM3Type , only : bgc_reaction_CENTURY_clm3_type
    use BGCReactionsCenturyCLMOType , only : bgc_reaction_CENTURY_clmo_type
    use abortutils                  , only : endrun
    use clm_varctl                  , only : iulog
    use tracer_varcon               , only : is_active_betr_bgc, do_betr_leaching

    ! !ARGUMENTS:
    class(bgc_reaction_type), allocatable :: bgc_reaction
    character(len=*), intent(in)          :: method
    character(len=*), parameter           :: subname = 'ctreate_bgc_reaction_type'

    select case(trim(method))
    case ("mock_run")
       allocate(bgc_reaction, source=bgc_reaction_mock_run_type())
    case ("century_bgc")
       is_active_betr_bgc = .true.
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_type())
    case ("century_bgcclm")
       is_active_betr_bgc = .true.
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_clm_type())
    case ("century_bgcECA")
       is_active_betr_bgc = .true.
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_ECA_type())
    case ("century_bgcclm3")
       is_active_betr_bgc=.true.
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_clm3_type())
    case ("century_bgcclmo")
       is_active_betr_bgc=.true.
       allocate(bgc_reaction, source=bgc_reaction_CENTURY_clmo_type())
    case ("betr_sminn")
       !this must be used together with clm45bgc
       do_betr_leaching = .true.
       allocate(bgc_reaction, source=bgc_reaction_sminn_type())
       !case ("o18_istope")    ! on hold
       !  allocate(bgc_reaction, source=bgc_reaction_O18ISO_type())
    case default
       write(iulog,*)subname //' ERROR: unknown method: ', method
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
  end function ctreate_bgc_reaction_type

end module BGCReactionsFactoryMod
