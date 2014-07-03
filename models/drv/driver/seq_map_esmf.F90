module seq_map_esmf
#ifdef USE_ESMF_LIB

! !USES:
  use esmf
  use esmfshr_mod
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: seq_map_esmf_init
  public :: seq_map_esmf_run
  public :: seq_map_esmf_final
  public :: seq_map_esmf_register

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!
! Author: Fei Liu
! ESMF compliant coupler  component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine seq_map_esmf_register(comp, rc)

    implicit none

    type(ESMF_GridComp), intent(inout)  :: comp
    integer,             intent(out)    :: rc

    rc = ESMF_SUCCESS

    print *, "In seq map register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      seq_map_esmf_init, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      seq_map_esmf_run, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      seq_map_esmf_final, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_map_esmf_init
!
! !DESCRIPTION:
!     initialize dead seq_map model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_map_esmf_init(comp, import_state, export_state, EClock, rc)

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

!EOP

    !--- local variables ---

    type(ESMF_Array)             :: gindex_s, gindex_d, factorArray, factorIndexArray
    type(ESMF_Array)             :: array_s, array_d
    type(ESMF_Routehandle)       :: routehandle
    integer, pointer             :: gindex(:), factorIndexList(:,:)
    real(ESMF_KIND_R8), pointer  :: factorList(:)
    character(len=64)            :: string = 'seq_map_esmf_init'

    !----------------------------
    ! Initial Setup
    !----------------------------
    ! Get temporary Arrays from States
    call ESMF_StateGet(import_state, itemName='gindex_s', array=gindex_s, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName='gindex_d', array=gindex_d, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName='factorArray', array=factorArray, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName='factorIndexArray', array=factorIndexArray, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Using gindex lists to set up Arrays for routehandle calculation
    ! These temporary Arrays are 2D of the size (1, size(gindex))
    call ESMF_ArrayGet(gindex_s, farrayPtr=gindex, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    array_s = mct2esmf_init(gindex, 'temparray_s', name='temparray_s', rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(gindex_d, farrayPtr=gindex, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    array_d = mct2esmf_init(gindex, 'temparray_d', name='temparray_d', rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(factorArray, farrayPtr=factorList, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(factorIndexArray, farrayPtr=factorIndexList, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! SMM Store
    call ESMF_ArraySMMStore(srcArray=array_s, dstArray=array_d, routehandle=routehandle, &
          factorList=factorList, factorIndexList=factorIndexList, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Attach routehandle to export State
    call ESMF_RoutehandleSet(routehandle, name="routehandle", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateAdd(export_state, (/routehandle/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Attach temporary Arrays to import/export State
    call ESMF_StateAdd(import_state, (/array_s/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateAdd(export_state, (/array_d/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

end subroutine seq_map_esmf_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_map_esmf_run
!
! !DESCRIPTION:
!     run method for dead seq_map model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_map_esmf_run(comp, import_state, export_state, EClock, rc)

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

!EOP

    !--- local variables ---

    type(ESMF_Array)             :: array_s, array_d
    type(ESMF_Routehandle)       :: routehandle

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Retrieve src and dst Arrays and routehandles
    call ESMF_StateGet(import_state, itemName='array_s', array=array_s, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName='array_d', array=array_d, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName='routehandle', routehandle=routehandle, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! SMM execution
    call ESMF_ArraySMM(srcArray=array_s, dstArray=array_d, routehandle=routehandle, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine seq_map_esmf_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_map_esmf_final
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_map_esmf_final(comp, import_state, export_state, EClock, rc)

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
!EOP

end subroutine seq_map_esmf_final
!===============================================================================

#endif
end module seq_map_esmf
