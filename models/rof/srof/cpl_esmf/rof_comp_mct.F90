module rof_comp_mct

! !USES:
  use shr_kind_mod, only: IN=>shr_kind_IN
  use seq_infodata_mod
  use seq_timemgr_mod
  use mct_mod
  use ESMF

  use seq_cdata_mod
  use seq_comm_mct     , only: seq_comm_getcompstates

  use esmfshr_mod
  use rof_comp_esmf
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: rof_init_mct
  public :: rof_run_mct
  public :: rof_final_mct
  public :: rof_register

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine rof_register(rof_petlist, ccsmComp, rof_comp, import_state, export_state)

  implicit none

  ! Arguments
  integer, pointer                  :: rof_petlist(:)
  type(ESMF_CplComp)                :: ccsmComp
  type(ESMF_GridComp),intent(out)   :: rof_comp
  type(ESMF_State)   ,intent(out)   :: import_state, export_state

  ! Local variables
  integer            :: rc
  
  ! Create Gridded Component and import/export States
  rof_comp = ESMF_GridCompCreate(name="rof_comp", petList=rof_petlist, rc=rc)
  if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create rof comp')
  call ESMF_GridCompSetServices(rof_comp, rof_register_esmf, rc=rc)
  if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register rof comp')
  import_state = ESMF_StateCreate(name="rof import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
  if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import rof state')
  export_state = ESMF_StateCreate(name="rof export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
  if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export rof state')
  
  ! Link attribute
  call ESMF_AttributeLink(ccsmComp, rof_comp, rc=rc)
  if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')
  
end subroutine

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_init_mct
!
! !DESCRIPTION:
!     initialize stub rof model
!
! !REVISION HISTORY:
!
! !Authors:
!   Fei Liu
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine rof_init_mct( EClock, cdata, x2d, d2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
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
    type(ESMF_GridComp)                   :: rof_comp
!-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)
    call seq_comm_getcompstates(COMPID, rof_comp, import_state, export_state)

    ! Copy infodata to state
    call esmfshr_infodata_infodata2state(infodata,export_state,ID=COMPID,rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! call into ESMF init method which calls ESMF run method
    call ESMF_GridCompInitialize(rof_comp, importState=import_state, exportState=export_state, &
         clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! copy export_state to infodata
    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine rof_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_run_mct
!
! !DESCRIPTION:
!     run method for stub rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine rof_run_mct(EClock, cdata, x2d, d2x)

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> stub
    type(mct_aVect)             ,intent(inout) :: d2x        ! stub   -> driver

!EOP
    type(seq_infodata_type)     , pointer :: infodata
    integer                               :: rc, urc, COMPID
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: rof_comp

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, infodata=infodata, ID=COMPID)
    call seq_comm_getcompstates(COMPID, rof_comp, import_state, export_state)

    call ESMF_AttributeSet(export_state, name="ID", value=COMPID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompRun(rof_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine rof_run_mct

!===============================================================================

subroutine rof_final_mct( EClock, cdata, x2d, d2x)

    implicit none

    !----- arguments -----

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d
    type(mct_aVect)             ,intent(inout) :: d2x

    !----- local -----
    integer             :: rc, urc, COMPID
    type(ESMF_State)    :: import_state, export_state
    type(ESMF_GridComp) :: rof_comp
    !----------------------------------------------------------------------------

    call seq_cdata_setptrs(cdata, ID=COMPID)
    call seq_comm_getcompstates(COMPID, rof_comp, import_state, export_state)

    !----------------------------------------------------------------------------
    ! Finalize routine 
    !----------------------------------------------------------------------------

    call ESMF_GridCompFinalize(rof_comp, importState=import_state, exportState=export_state, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! destroy component and states
    call ESMF_StateDestroy(import_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateDestroy(export_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompDestroy(rof_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine rof_final_mct
!===============================================================================

end module rof_comp_mct
