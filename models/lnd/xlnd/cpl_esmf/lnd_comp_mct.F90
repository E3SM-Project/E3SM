module lnd_comp_mct

! !USES:
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS
  use seq_infodata_mod
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

!
! Author: Fei Liu
! This module is a wrapper layer between ccsm driver and ESMF dead lnd component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine lnd_register(lnd_petlist, ccsmComp, lnd_comp, import_state, export_state)

    implicit none

    integer, pointer                  :: lnd_petlist(:)
    type(ESMF_CplComp) ,intent(inout) :: ccsmComp
    type(ESMF_GridComp),intent(out)   :: lnd_comp
    type(ESMF_State)   ,intent(out)   :: import_state, export_state

    integer            :: rc

    lnd_comp = ESMF_GridCompCreate(name="lnd_comp", petList=lnd_petlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create lnd comp')
    call ESMF_GridCompSetServices(lnd_comp, lnd_register_esmf, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register lnd comp')
    import_state = ESMF_StateCreate(name="lnd import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import lnd state')
    export_state = ESMF_StateCreate(name="lnd export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export lnd state')

    call ESMF_AttributeLink(ccsmComp, lnd_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

end subroutine

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_init_mct
!
! !DESCRIPTION:
!     initialize dead lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_init_mct( EClock, cdata, x2d, d2x, &
                               cdata_s, x2s, s2x, NLFilename )

    implicit none

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
    type(mct_gsMap)             , pointer :: gsMap_s
    type(mct_gGrid)             , pointer :: dom_s
    type(seq_infodata_type)     , pointer :: infodata
    type(ESMF_Array)                      :: d2x_a, x2d_a, dom_a
    type(ESMF_Array)                      :: s2x_a, x2s_a, doms_a
    integer                               :: phase, rc, urc
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: lnd_comp
!-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)
    call seq_comm_getcompstates(COMPID, lnd_comp, import_state, export_state)

    call seq_cdata_setptrs(cdata_s, &
         gsMap=gsMap_s, dom=dom_s)

    ! We are running into a coupling issue here...things are updated externally
    ! in infodata but not states; Init during phase 1 only.
    call seq_infodata_getData(infodata, lnd_phase=phase)
    if (phase > 1) then
        return
    endif

    ! Copy infodata to state

    call esmfshr_infodata_infodata2state(infodata,export_state,ID=COMPID,rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! call into ESMF init method
    call ESMF_GridCompInitialize(lnd_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! copy export_state to infodata

    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! extract d2x and x2d esmf arrays
    call ESMF_StateGet(export_state, itemName="domain", array=dom_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(export_state, itemName="domain_s", array=doms_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(export_state, itemName="d2x", array=d2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(export_state, itemName="s2x", array=s2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(import_state, itemName="x2d", array=x2d_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(import_state, itemName="x2s", array=x2s_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize domains

    call esmf2mct_init(d2x_a, COMPID, gsMap, mpicom, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_init(s2x_a, COMPID, gsMap_s, mpicom, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Set domain and state arrays

    call esmf2mct_init(dom_a, dom, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(dom_a, dom%data, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(doms_a, dom_s, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(doms_a, dom_s%data, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(x2d_a, x2d, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call mct_aVect_zero(x2d)

    call esmf2mct_init(x2s_a, x2s, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call mct_aVect_zero(x2s)

    call esmf2mct_init(d2x_a, d2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(d2x_a, d2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(s2x_a, s2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(s2x_a, s2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine lnd_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_run_mct
!
! !DESCRIPTION:
!     run method for dead lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_run_mct( EClock, cdata,  x2d, d2x, cdata_s, x2s, s2x)

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
    integer(IN)                           :: COMPID
    type(seq_infodata_type)     , pointer :: infodata
    type(ESMF_Array)                      :: d2x_a, s2x_a
    integer                               :: rc, urc
    logical  :: log1, log2
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: lnd_comp

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID, infodata=infodata)
    call seq_comm_getcompstates(COMPID, lnd_comp, import_state, export_state)

    ! update lnd state according to infodata from coupler
    call esmfshr_infodata_infodata2state(infodata, export_state, ID=COMPID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompRun(lnd_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! convert state back to infodata
    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! copy values from mct to esmf

    call ESMF_StateGet(export_state, itemName="d2x", array=d2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(d2x_a, d2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName="s2x", array=s2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(s2x_a, s2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
end subroutine lnd_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_final_mct
!
! !DESCRIPTION:
!     finalize method for dead model
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
