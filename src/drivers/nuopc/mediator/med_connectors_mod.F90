module med_connectors_mod

  !-----------------------------------------------------------------------------
  ! Connector phases 
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use esmFlds               , only: compatm, compocn, compice
  use esmFlds               , only: complnd, comprof, compwav, compglc
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_State_diagnose
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_State_reset
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_copy
  use med_infodata_mod      , only: med_infodata_CopyStateToInfodata
  use med_infodata_mod      , only: med_infodata_CopyInfodataToState
  use med_infodata_mod      , only: med_infodata
  use med_internalstate_mod
  use med_constants_mod

  implicit none
  private

  integer                       :: dbug_flag = med_constants_dbug_flag
  logical                       :: statewrite_flag = med_constants_statewrite_flag
  real(ESMF_KIND_R8), parameter :: spval = med_constants_spval
  real(ESMF_KIND_R8), parameter :: czero = med_constants_czero
  character(*)      , parameter :: u_FILE_u = __FILE__
  integer                       :: dbrc

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_connectors_prep_med2atm
  public med_connectors_prep_med2ocn
  public med_connectors_prep_med2ice
  public med_connectors_prep_med2lnd
  public med_connectors_prep_med2rof
  public med_connectors_prep_med2wav
  public med_connectors_prep_med2glc
  public med_connectors_post_atm2med
  public med_connectors_post_ocn2med
  public med_connectors_post_ice2med
  public med_connectors_post_lnd2med
  public med_connectors_post_rof2med
  public med_connectors_post_wav2med
  public med_connectors_post_glc2med

  !--------------------------------------------------------------------------
  ! Private
  !--------------------------------------------------------------------------

  private med_connectors_prep_generic
  private med_connectors_post_generic
  private med_connectors_diagnose

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_connectors_prep_generic(gcomp, type, rc)
    type(ESMF_GridComp)  :: gcomp
    character(len=*), intent(in) :: type
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(InternalState)        :: is_local
    logical                    :: diagnose
    logical                    :: connected
    integer                    :: n
    character(len=*),parameter :: subname='(med_connectors_prep_generic)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(type)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return


    !-------------------------
    ! diagnose export state
    ! update scalar data in Exp and Imp State
    !-------------------------

    select case (type)

    case('atm')
      is_local%wrap%conn_prep_cnt(compatm) = is_local%wrap%conn_prep_cnt(compatm) + 1
      call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(compatm), value=spval, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(compatm), is_local%wrap%FBExp(compatm), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_connectors_diagnose(is_local%wrap%NStateExp(compatm), is_local%wrap%conn_prep_cnt(compatm), "med_to_atm", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(compatm),'cpl2atm',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(compatm),'cpl2atm',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('ocn')
      is_local%wrap%conn_prep_cnt(compocn) = is_local%wrap%conn_prep_cnt(compocn) + 1
      call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(compocn), value=spval, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(compocn), is_local%wrap%FBExp(compocn), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_connectors_diagnose(is_local%wrap%NStateExp(compocn), is_local%wrap%conn_prep_cnt(compocn), "med_to_ocn", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(compocn),'cpl2ocn',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(compocn),'cpl2ocn',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('ice')
      is_local%wrap%conn_prep_cnt(compice) = is_local%wrap%conn_prep_cnt(compice) + 1
      call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(compice), value=spval, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(compice), is_local%wrap%FBExp(compice), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_connectors_diagnose(is_local%wrap%NStateExp(compice), is_local%wrap%conn_prep_cnt(compice), "med_to_ice", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(compice),'cpl2ice',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(compice),'cpl2ice',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('lnd')
      is_local%wrap%conn_prep_cnt(complnd) = is_local%wrap%conn_prep_cnt(complnd) + 1
      call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(complnd), value=spval, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(complnd), is_local%wrap%FBExp(complnd), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_connectors_diagnose(is_local%wrap%NStateExp(complnd), is_local%wrap%conn_prep_cnt(complnd), "med_to_lnd", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(complnd),'cpl2lnd',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(complnd),'cpl2lnd',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('rof')
      is_local%wrap%conn_prep_cnt(comprof) = is_local%wrap%conn_prep_cnt(comprof) + 1
      call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(comprof), value=spval, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(comprof), is_local%wrap%FBExp(comprof), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_connectors_diagnose(is_local%wrap%NStateExp(comprof), is_local%wrap%conn_prep_cnt(comprof), "med_to_rof", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(comprof),'cpl2rof',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(comprof),'cpl2rof',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('wav')
      is_local%wrap%conn_prep_cnt(compwav) = is_local%wrap%conn_prep_cnt(compwav) + 1
      call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(compwav), value=spval, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(compwav), is_local%wrap%FBExp(compwav), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_connectors_diagnose(is_local%wrap%NStateExp(compwav), is_local%wrap%conn_prep_cnt(compwav), "med_to_wav", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(compwav),'cpl2wav',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(compwav),'cpl2wav',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('glc')
      is_local%wrap%conn_prep_cnt(compglc) = is_local%wrap%conn_prep_cnt(compglc) + 1
      call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(compglc), value=spval, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(compglc), is_local%wrap%FBExp(compglc), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_connectors_diagnose(is_local%wrap%NStateExp(compglc), is_local%wrap%conn_prep_cnt(compglc), "med_to_glc", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(compglc),'cpl2glc',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(compglc),'cpl2glc',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case default
      rc = ESMF_Failure
      call ESMF_LogWrite(trim(subname)//trim(type)//" unsupported", ESMF_LOGMSG_INFO, rc=dbrc)

    end select

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(type)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_generic(gcomp, type, rc)
    type(ESMF_GridComp)           :: gcomp
    character(len=*), intent(in)  :: type
    integer,          intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(InternalState)        :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_generic)'

    ! Note: for information obtained by the mediator always write out the state
    ! if statewrite_flag is .true.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(type)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------
    ! diagnose import state
    ! copy import state scalar data to local datatype
    !-------------------------

    select case (type)

    case('atm')
      is_local%wrap%conn_post_cnt(compatm) = is_local%wrap%conn_post_cnt(compatm) + 1
      call med_connectors_diagnose(is_local%wrap%NStateImp(compatm), is_local%wrap%conn_post_cnt(compatm), " med_from_atm", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(compatm),med_infodata,'atm2cpl',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compatm,compatm), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(compatm,compatm), is_local%wrap%NStateImp(compatm), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('ocn')
      is_local%wrap%conn_post_cnt(compocn) = is_local%wrap%conn_post_cnt(compocn) + 1
      call med_connectors_diagnose(is_local%wrap%NStateImp(compocn), is_local%wrap%conn_post_cnt(compocn), " med_from_ocn", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(compocn),med_infodata,'ocn2cpl',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compocn,compocn), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(compocn,compocn), is_local%wrap%NStateImp(compocn), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('ice')
      is_local%wrap%conn_post_cnt(compice) = is_local%wrap%conn_post_cnt(compice) + 1
      call med_connectors_diagnose(is_local%wrap%NStateImp(compice), is_local%wrap%conn_post_cnt(compice), " med_from_ice", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(compice),med_infodata,'ice2cpl',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compice,compice), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(compice,compice), is_local%wrap%NStateImp(compice), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('lnd')
      is_local%wrap%conn_post_cnt(complnd) = is_local%wrap%conn_post_cnt(complnd) + 1
      call med_connectors_diagnose(is_local%wrap%NStateImp(complnd), is_local%wrap%conn_post_cnt(complnd), " med_from_lnd", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(complnd),med_infodata,'lnd2cpl',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(complnd,complnd), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(complnd,complnd), is_local%wrap%NStateImp(complnd), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('rof')
      is_local%wrap%conn_post_cnt(comprof) = is_local%wrap%conn_post_cnt(comprof) + 1
      call med_connectors_diagnose(is_local%wrap%NStateImp(comprof), is_local%wrap%conn_post_cnt(comprof), " med_from_rof", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(comprof),med_infodata,'rof2cpl',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(comprof,comprof), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(comprof,comprof), is_local%wrap%NStateImp(comprof), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('wav')
      is_local%wrap%conn_post_cnt(compwav) = is_local%wrap%conn_post_cnt(compwav) + 1
      call med_connectors_diagnose(is_local%wrap%NStateImp(compwav), is_local%wrap%conn_post_cnt(compwav), " med_from_wav", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(compwav),med_infodata,'wav2cpl',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compwav,compwav), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(compwav,compwav), is_local%wrap%NStateImp(compwav), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case('glc')
      is_local%wrap%conn_post_cnt(compglc) = is_local%wrap%conn_post_cnt(compglc) + 1
      call med_connectors_diagnose(is_local%wrap%NStateImp(compglc), is_local%wrap%conn_post_cnt(compglc), " med_from_glc", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(compglc),med_infodata,'glc2cpl',is_local%wrap%mpicom,rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compglc,compglc), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(compglc,compglc), is_local%wrap%NStateImp(compglc), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case default
      rc = ESMF_Failure
      call ESMF_LogWrite(trim(subname)//trim(type)//" unsupported", ESMF_LOGMSG_INFO, rc=dbrc)

    end select

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(type)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2atm(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_prep_med2atm)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'atm', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_med2atm

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ocn(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_prep_med2ocn)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'ocn', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_med2ocn

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ice(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_prep_med2ice)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'ice', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_med2ice

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2lnd(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_prep_med2lnd)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'lnd', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_med2lnd

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2rof(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_prep_med2rof)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'rof', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_med2rof

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2wav(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_prep_med2wav)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'wav', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_med2wav

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2glc(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_prep_med2glc)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'glc', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_prep_med2glc

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_atm2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_atm2med)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'atm', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_atm2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ocn2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_ocn2med)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'ocn', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_ocn2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ice2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_ice2med)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'ice', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_ice2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_lnd2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_lnd2med)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'lnd', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_lnd2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_rof2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_rof2med)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'rof', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_rof2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_wav2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_wav2med)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'wav', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_wav2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_glc2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    character(len=*),parameter :: subname='(med_connectors_post_glc2med)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'glc', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_glc2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_diagnose(State, cntr, string, rc)
    type(ESMF_State), intent(in)    :: State
    integer         , intent(inout) :: cntr
    character(len=*), intent(in)    :: string
    integer         , intent(out)   :: rc

    ! local variables
    integer :: fieldCount
    character(ESMF_MAXSTR),pointer :: fieldnamelist(:)
    character(len=*),parameter :: subname='(med_connectors_diagnose)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(string)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Obtain the field names in State - allocate memory which will be deallocated at the end
    allocate(fieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=fieldnamelist, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
      call shr_nuopc_methods_State_diagnose(State, string=trim(subname)//trim(string), rc=rc)
    endif

    ! Write out the fields in State to netcdf files
    if (cntr > 0 .and. statewrite_flag) then
      call ESMF_LogWrite(trim(subname)//trim(string)//": writing out fields", ESMF_LOGMSG_INFO, rc=dbrc)
      call NUOPC_Write(State, &
        fieldnamelist(1:fieldCount), &
        "field_"//trim(string)//"_", timeslice=cntr, &
        overwrite=.true., relaxedFlag=.true., rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    deallocate(fieldnamelist)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(string)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_diagnose

  !-----------------------------------------------------------------------------

end module med_connectors_mod
