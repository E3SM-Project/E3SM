module glc_comp_mct

! !USES:
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS
  use seq_infodata_mod
  use mct_mod
  use ESMF

  use seq_cdata_mod
  use seq_comm_mct     , only: seq_comm_getcompstates

  use esmfshr_mod
  use glc_comp_esmf

! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: glc_init_mct
  public :: glc_run_mct
  public :: glc_final_mct
  public :: glc_register

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

   type(ESMF_GridComp)     :: glc_comp
   type(ESMF_State)        :: import_state, export_state
!
! Author: Fei Liu
! This module is a wrapper layer between ccsm driver and ESMF dead glc component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine glc_register(glc_petlist, ccsmComp, glc_comp, import_state, export_state)

    implicit none

    integer, pointer                  :: glc_petlist(:)
    type(ESMF_CplComp) ,intent(inout) :: ccsmComp
    type(ESMF_GridComp),intent(out)   :: glc_comp
    type(ESMF_State)   ,intent(out)   :: import_state, export_state

    integer            :: rc

    glc_comp = ESMF_GridCompCreate(name="glc_comp", petList=glc_petlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create glc comp')
    call ESMF_GridCompSetServices(glc_comp, glc_register_esmf, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register glc comp')
    import_state = ESMF_StateCreate(name="glc import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import glc state')
    export_state = ESMF_StateCreate(name="glc export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export glc state')

    call ESMF_AttributeLink(ccsmComp, glc_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

end subroutine

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_init_mct
!
! !DESCRIPTION:
!     initialize dead glc model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine glc_init_mct( EClock, cdata, x2d, d2x, NLFilename )

    implicit none

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
    type(ESMF_Array)                      :: d2x_a, x2d_a, dom_a
    integer                               :: phase, rc, urc
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: glc_comp
!-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)

    call seq_comm_getcompstates(COMPID, glc_comp, import_state, export_state)

    ! We are running into a coupling issue here...things are updated externally
    ! in infodata but not states; Init during phase 1 only.
    call seq_infodata_getData(infodata, glc_phase=phase)
    if (phase > 1) then
        call ESMF_AttributeSet(export_state, name="glc_phase", value=phase, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
        return
    endif

    ! Copy infodata to state

    call esmfshr_infodata_infodata2state(infodata,export_state,ID=COMPID,rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! call into ESMF init method
    call ESMF_GridCompInitialize(glc_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! copy export_state to infodata
    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! extract d2x and x2d esmf arrays
    call ESMF_StateGet(export_state, itemName="domain", array=dom_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(export_state, itemName="d2x", array=d2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(import_state, itemName="x2d", array=x2d_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize MCT domain
    call esmf2mct_init(d2x_a, COMPID, gsMap, mpicom, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(dom_a, dom, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(dom_a, dom%data, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize MCT attribute vectors
    call esmf2mct_init(x2d_a, x2d, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call mct_aVect_zero(x2d)

    call esmf2mct_init(d2x_a, d2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    ! copy values back to d2x
    call esmf2mct_copy(d2x_a, d2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine glc_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_run_mct
!
! !DESCRIPTION:
!     run method for dead glc model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine glc_run_mct(EClock, cdata, x2d, d2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

!EOP
    integer(IN)                           :: COMPID
    type(seq_infodata_type)     , pointer :: infodata
    type(ESMF_Array)                      :: d2x_a
    integer                               :: rc, urc
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: glc_comp

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID, infodata=infodata)

    call seq_comm_getcompstates(COMPID, glc_comp, import_state, export_state)

    call esmfshr_infodata_infodata2state(infodata, export_state, ID=COMPID,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompRun(glc_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! convert state back to infodata
    call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! copy values back to d2x
    call ESMF_StateGet(export_state, itemName="d2x", array=d2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(d2x_a, d2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
end subroutine glc_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_final_mct
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine glc_final_mct(EClock, cdata, x2d, d2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

!EOP

    integer                               :: rc, urc, COMPID
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: glc_comp

    !----------------------------------------------------------------------------
    ! Finalize routine 
    !----------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID)

    call seq_comm_getcompstates(COMPID, glc_comp, import_state, export_state)

    call ESMF_GridCompFinalize(glc_comp, importState=import_state, exportState=export_state, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! destroy component and states
    call ESMF_StateDestroy(import_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateDestroy(export_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompDestroy(glc_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine glc_final_mct
!===============================================================================

end module glc_comp_mct
