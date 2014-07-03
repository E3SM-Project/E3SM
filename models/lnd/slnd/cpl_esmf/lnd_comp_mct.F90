module lnd_comp_mct

! !USES:
  use shr_kind_mod, only: IN=>shr_kind_IN
  use seq_infodata_mod
  use seq_timemgr_mod
  use mct_mod
  use ESMF

  use seq_cdata_mod
  use seq_comm_mct     , only: seq_comm_getcompstates

  use esmfshr_mod
  use lnd_comp_esmf
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct
  public :: lnd_register

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine lnd_register(lnd_petlist, ccsmComp, lnd_comp, import_state, export_state)

   implicit none

   ! Arguments
   integer, pointer                  :: lnd_petlist(:)
   type(ESMF_CplComp)                :: ccsmComp
   type(ESMF_GridComp),intent(out)   :: lnd_comp
   type(ESMF_State)   ,intent(out)   :: import_state, export_state

  ! Local variables
   integer            :: rc

  ! Create Gridded Component and import/export States
   lnd_comp = ESMF_GridCompCreate(name="lnd_comp", petList=lnd_petlist, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create lnd comp')
   call ESMF_GridCompSetServices(lnd_comp, lnd_register_esmf, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register lnd comp')
   import_state = ESMF_StateCreate(name="lnd import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import lnd state')
   export_state = ESMF_StateCreate(name="lnd export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export lnd state')

   ! Link attribute
   call ESMF_AttributeLink(ccsmComp, lnd_comp, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

end subroutine

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_init_mct
!
! !DESCRIPTION:
!     initialize stub lnd model
!
! !REVISION HISTORY:
!
! !Authors:
!   Fei Liu
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine lnd_init_mct( EClock, cdata, x2d, d2x, cdata_s, x2s, s2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    type(seq_cdata)             , intent(inout) :: cdata_s
    type(mct_aVect)             , intent(inout) :: x2s, s2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    !--- local variables ---
    integer(IN)                           :: mpicom
    integer(IN)                           :: COMPID
    type(mct_gsMap)             , pointer :: gsMap
    type(mct_gGrid)             , pointer :: dom
    type(seq_infodata_type)     , pointer :: infodata
    integer                               :: rc, urc
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: lnd_comp
!-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)
    call seq_comm_getcompstates(COMPID, lnd_comp, import_state, export_state)

    ! Copy infodata to state
    call esmfshr_infodata_infodata2state(infodata,export_state,ID=COMPID,rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! call into ESMF init method which calls ESMF run method
    call ESMF_GridCompInitialize(lnd_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! copy export_state to infodata
    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine lnd_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_run_mct
!
! !DESCRIPTION:
!     run method for stub lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_run_mct(EClock, cdata, x2d, d2x, cdata_s, x2s, s2x)

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> stub
    type(mct_aVect)             ,intent(inout) :: d2x        ! stub   -> driver
    type(seq_cdata)             ,intent(inout) :: cdata_s
    type(mct_aVect)             ,intent(inout) :: x2s        ! driver -> stub
    type(mct_aVect)             ,intent(inout) :: s2x        ! stub   -> driver

!EOP
    type(seq_infodata_type)     , pointer :: infodata
    integer                               :: rc, urc, COMPID
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: lnd_comp

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, infodata=infodata, ID=COMPID)
    call seq_comm_getcompstates(COMPID, lnd_comp, import_state, export_state)

    call ESMF_AttributeSet(export_state, name="ID", value=COMPID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompRun(lnd_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine lnd_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_final_mct
!
! !DESCRIPTION:
!     finalize method for stub model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_final_mct( EClock, cdata, x2d, d2x, cdata_s, x2s, s2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver
    type(seq_cdata)             ,intent(inout) :: cdata_s
    type(mct_aVect)             ,intent(inout) :: x2s
    type(mct_aVect)             ,intent(inout) :: s2x

!EOP

    integer             :: rc, urc, COMPID
    type(ESMF_State)    :: import_state, export_state
    type(ESMF_GridComp) :: lnd_comp

    !----------------------------------------------------------------------------
    ! Finalize routine 
    !----------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID)
    call seq_comm_getcompstates(COMPID, lnd_comp, import_state, export_state)

    call ESMF_GridCompFinalize(lnd_comp, importState=import_state, exportState=export_state, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! destroy component and states
    call ESMF_StateDestroy(import_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateDestroy(export_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompDestroy(lnd_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine lnd_final_mct
!===============================================================================

end module lnd_comp_mct
