module med_connectors_mod

  !-----------------------------------------------------------------------------
  ! Connector phases
  !-----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_Failure
  use ESMF                  , only : ESMF_State, ESMF_Clock, ESMF_GridComp
  use med_internalstate_mod , only : InternalState
  use shr_nuopc_utils_mod , only : shr_nuopc_utils_ChkErr
  use med_constants_mod     , only : spval => med_constants_spval
  use med_constants_mod     , only : czero => med_constants_czero

  implicit none
  private
  character(*)      , parameter :: u_FILE_u = &
       __FILE__

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

  subroutine med_connectors_prep_generic(gcomp, type, compid, rc)
    use ESMF                  , only : ESMF_GridCompGet, ESMF_VMGet
    use med_infodata_mod      , only : med_infodata_CopyStateToInfodata
    use med_infodata_mod      , only : med_infodata_CopyInfodataToState
    use med_infodata_mod      , only : med_infodata
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_copy
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)          :: gcomp
    character(len=*), intent(in) :: type
    integer, intent(in)          :: compid
    integer, intent(out)         :: rc

    ! local variables
    type(ESMF_Clock)    :: clock
    type(InternalState) :: is_local
    logical             :: diagnose
    logical             :: connected
    integer             :: n
    integer             :: dbrc
    integer             :: mytask
    character(len=10)   :: med2comp
    character(len=7)    :: cpl2comp
    character(len=*),parameter :: subname='(med_connectors_prep_generic)'
    !---------------------------------------------
    call t_startf('MED:'//subname)
    call ESMF_LogWrite(trim(subname)//trim(type)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(is_local%wrap%vm, localPet=mytask, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    !-------------------------
    ! diagnose export state
    ! update scalar data in Exp and Imp State
    !-------------------------
    med2comp = "med_to_"//type
    cpl2comp = "cpl2"//type

    is_local%wrap%conn_prep_cnt(compid) = is_local%wrap%conn_prep_cnt(compid) + 1
    call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(compid), value=spval, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(compid), is_local%wrap%FBExp(compid), rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call med_connectors_diagnose(is_local%wrap%NStateExp(compid), is_local%wrap%conn_prep_cnt(compid), med2comp, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(compid), cpl2comp, mytask, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(compid), cpl2comp, mytask, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//trim(type)//": done", ESMF_LOGMSG_INFO, rc=rc)

    call t_stopf('MED:'//subname)

  end subroutine med_connectors_prep_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_generic(gcomp, type, compid, rc)

    use ESMF                  , only : ESMF_GridCompGet
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_copy
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use med_infodata_mod      , only : med_infodata
    use med_infodata_mod      , only : med_infodata_CopyStateToInfodata
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)           :: gcomp
    character(len=*), intent(in)  :: type
    integer,          intent(in)  :: compid
    integer,          intent(out) :: rc

    ! local variables
    type(ESMF_Clock)    :: clock
    type(InternalState) :: is_local
    integer             :: dbrc
    character(len=10)   :: comp2med
    character(len=7)    :: comp2cpl
    character(len=*),parameter :: subname='(med_connectors_post_generic)'
    !---------------------------------------------

    ! Note: for information obtained by the mediator always write out the state
    ! if statewrite_flag is .true.
    rc = ESMF_SUCCESS
    call t_startf('MED:'//subname)

    call ESMF_LogWrite(trim(subname)//trim(type)//": called", ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------
    ! diagnose import state
    ! copy import state scalar data to local datatype
    !-------------------------
    comp2med = "med_from_"//type
    comp2cpl = type//"2cpl"

    is_local%wrap%conn_post_cnt(compid) = is_local%wrap%conn_post_cnt(compid) + 1
    call med_connectors_diagnose(is_local%wrap%NStateImp(compid), is_local%wrap%conn_post_cnt(compid),comp2med, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(compid),med_infodata, comp2cpl ,is_local%wrap%vm,rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compid,compid), value=czero, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(compid,compid), is_local%wrap%NStateImp(compid), rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//trim(type)//": done", ESMF_LOGMSG_INFO, rc=rc)

    call t_stopf('MED:'//subname)

  end subroutine med_connectors_post_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2atm(gcomp, rc)
    use perf_mod, only : t_startf, t_stopf
    use esmFlds,  only : compatm
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2atm)'
    !---------------------------------------------
    call t_startf('MED:'//subname)

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)

    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'atm', compatm, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)
    call t_stopf('MED:'//subname)

  end subroutine med_connectors_prep_med2atm

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ocn(gcomp, rc)
    use esmFlds,  only : compocn
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2ocn)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)

    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'ocn', compocn, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_prep_med2ocn

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ice(gcomp, rc)
    use esmFlds,  only : compice
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2ice)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'ice', compice, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_prep_med2ice

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2lnd(gcomp, rc)
    use esmFlds,  only : complnd
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2lnd)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'lnd', complnd, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_prep_med2lnd

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2rof(gcomp, rc)
    use esmFlds,  only : comprof
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer                       :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2rof)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'rof', comprof, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_prep_med2rof

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2wav(gcomp, rc)
    use esmFlds,  only : compwav
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer                       :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2wav)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)

    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'wav', compwav, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_prep_med2wav

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2glc(gcomp, rc)
    use esmFlds,  only : compglc
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2glc)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'glc', compglc, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_prep_med2glc

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_atm2med(gcomp, rc)
    use esmFlds,  only : compatm
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_atm2med)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)

    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'atm', compatm, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_post_atm2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ocn2med(gcomp, rc)
    use esmFlds,  only : compocn
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_ocn2med)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'ocn', compocn, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_post_ocn2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ice2med(gcomp, rc)
    use esmFlds,  only : compice
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_ice2med)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'ice', compice, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_post_ice2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_lnd2med(gcomp, rc)
    use esmFlds,  only : complnd
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_lnd2med)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'lnd', complnd, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_post_lnd2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_rof2med(gcomp, rc)
    use esmFlds,  only : comprof
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_rof2med)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'rof', comprof, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_post_rof2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_wav2med(gcomp, rc)
    use esmFlds,  only : compwav
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_wav2med)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)

    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'wav', compwav, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_post_wav2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_glc2med(gcomp, rc)
    use esmFlds,  only : compglc
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_glc2med)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, 'glc', compglc, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_post_glc2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_diagnose(State, cntr, string, rc)

    use ESMF                  , only : ESMF_State, ESMF_MAXSTR, ESMF_StateGet
    use NUOPC                 , only : NUOPC_Write
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_diagnose
    use med_constants_mod     , only : statewrite_flag => med_constants_statewrite_flag
    use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag

    ! input/output variables
    type(ESMF_State), intent(in)    :: State
    integer         , intent(inout) :: cntr
    character(len=*), intent(in)    :: string
    integer         , intent(out)   :: rc

    ! local variables
    integer                        :: fieldCount
    character(ESMF_MAXSTR),pointer :: fieldnamelist(:)
    integer                        :: dbrc
    character(len=*),parameter     :: subname='(med_connectors_diagnose)'
    !---------------------------------------------

    call ESMF_LogWrite(trim(subname)//trim(string)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! Obtain the field names in State - allocate memory which will be deallocated at the end
    allocate(fieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=fieldnamelist, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
      call shr_nuopc_methods_State_diagnose(State, string=trim(subname)//trim(string), rc=rc)
    endif

    ! Write out the fields in State to netcdf files
    if (cntr > 0 .and. statewrite_flag) then
      call ESMF_LogWrite(trim(subname)//trim(string)//": writing out fields", ESMF_LOGMSG_INFO, rc=rc)
      call NUOPC_Write(State, &
        fieldnamelist(1:fieldCount), &
        "field_"//trim(string)//"_", timeslice=cntr, &
        overwrite=.true., relaxedFlag=.true., rc=rc)
      if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    deallocate(fieldnamelist)

    call ESMF_LogWrite(trim(subname)//trim(string)//": done", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine med_connectors_diagnose

  !-----------------------------------------------------------------------------

end module med_connectors_mod
