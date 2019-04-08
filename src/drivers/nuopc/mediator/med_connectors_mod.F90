module med_connectors_mod

  !-----------------------------------------------------------------------------
  ! Connector phases
  !-----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF                  , only : ESMF_VM, ESMF_VMGet, ESMF_SUCCESS
  use med_internalstate_mod , only : InternalState
  use shr_nuopc_utils_mod   , only : shr_nuopc_utils_ChkErr
  use esmFlds               , only : compatm, compocn, compice, complnd, comprof, compwav, compglc
  use med_infodata_mod      , only : med_infodata 
  use med_infodata_mod      , only : med_infodata_CopyInfodataToState
  use med_infodata_mod      , only : med_infodata_CopyStateToInfodata

  implicit none
  private
  character(*), parameter :: u_FILE_u = &
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

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_connectors_prep_generic(gcomp, type, compid, rc)

    ! input/output variables
    type(ESMF_GridComp)          :: gcomp
    character(len=*), intent(in) :: type
    integer, intent(in)          :: compid
    integer, intent(out)         :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm 
    integer             :: mytask
    character(len=*),parameter :: subname='(med_connectors_prep_generic)'
    !---------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_infodata_CopyInfodataToState(med_infodata, is_local%wrap%NStateExp(compid), "cpl2"//type, mytask, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_generic(gcomp, type, compid, rc)

    ! input/output variables
    type(ESMF_GridComp)           :: gcomp
    character(len=*), intent(in)  :: type
    integer,          intent(in)  :: compid
    integer,          intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_VM)       :: vm
    character(len=*),parameter :: subname='(med_connectors_post_generic)'
    !---------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(compid), med_infodata, type//"2cpl", vm, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2atm(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_prep_med2atm)'
    !---------------------------------------------

    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, 'atm', compatm, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_med2atm

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ocn(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_prep_med2ocn)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_prep_generic(gcomp, 'ocn', compocn, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_med2ocn

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ice(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_prep_med2ice)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_prep_generic(gcomp, 'ice', compice, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_med2ice

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2lnd(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_prep_med2lnd)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_prep_generic(gcomp, 'lnd', complnd, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_med2lnd

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2rof(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_prep_med2rof)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_prep_generic(gcomp, 'rof', comprof, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_med2rof

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2wav(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_prep_med2wav)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_prep_generic(gcomp, 'wav', compwav, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_med2wav

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2glc(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_prep_med2glc)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_prep_generic(gcomp, 'glc', compglc, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_prep_med2glc

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_atm2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_post_atm2med)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_post_generic(gcomp, 'atm', compatm, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_atm2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ocn2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_post_ocn2med)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_post_generic(gcomp, 'ocn', compocn, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_ocn2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ice2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_post_ice2med)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_post_generic(gcomp, 'ice', compice, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_ice2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_lnd2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_post_lnd2med)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
   call med_connectors_post_generic(gcomp, 'lnd', complnd, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_lnd2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_rof2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_post_rof2med)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_post_generic(gcomp, 'rof', comprof, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_rof2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_wav2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_post_wav2med)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
   call med_connectors_post_generic(gcomp, 'wav', compwav, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_wav2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_glc2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter :: subname='(med_connectors_post_glc2med)'
    !---------------------------------------------

    rc = ESMF_SUCCESS
    call med_connectors_post_generic(gcomp, 'glc', compglc, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_connectors_post_glc2med

end module med_connectors_mod
