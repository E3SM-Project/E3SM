module rof_comp_mct

   use shr_kind_mod, only: CS=>shr_kind_CS, IN=>shr_kind_IN, R8=>shr_kind_R8

   use mct_mod
   use esmf
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_comm_mct     , only: seq_comm_getcompstates

   use esmfshr_mod
   use rof_comp_esmf

   implicit none

   public :: rof_init_mct
   public :: rof_run_mct
   public :: rof_final_mct
   public :: rof_register

   private ! except

   save ! save everything

!
! Author: Mariana Vertenstein
! This module is a wrapper layer between ccsm driver and ESMF data rof component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

subroutine rof_register(rof_petlist, ccsmComp, rof_comp, import_state, export_state)

   implicit none

   integer, pointer                  :: rof_petlist(:)
   type(ESMF_CplComp)                :: ccsmComp
   type(ESMF_GridComp),intent(out)   :: rof_comp
   type(ESMF_State), intent(out)     :: import_state
   type(ESMF_State), intent(out)     :: export_state

   integer            :: rc

   rof_comp = ESMF_GridCompCreate(name="rof_comp", petList=rof_petlist, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create lnd comp')

   call ESMF_GridCompSetServices(rof_comp, rof_register_esmf, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register lnd comp')

   import_state = ESMF_StateCreate(name="lnd import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import lnd state')

   export_state = ESMF_StateCreate(name="lnd export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export lnd state')

   call ESMF_AttributeLink(ccsmComp, rof_comp, rc=rc)
   if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

end subroutine

!===============================================================================

  subroutine rof_init_mct( EClock, cdata, x2d, d2x, NLFilename )

   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_Clock),intent(inout)              :: EClock
   type(seq_cdata), intent(inout)              :: cdata
   type(mct_aVect), intent(inout)              :: x2d   
   type(mct_aVect), intent(inout)              :: d2x   
   character(len=*), optional,   intent(in)    :: NLFilename ! Namelist filename
   
   !----- local -----
   integer(IN)                           :: COMPID
   integer(IN)                           :: mpicom
   type(mct_gsMap)        , pointer      :: gsmap
   type(mct_gGrid)        , pointer      :: dom
   type(seq_infodata_type), pointer      :: infodata
   integer                               :: rc, urc
   integer(IN)                           :: phase
   type(ESMF_State)                      :: import_state, export_state
   type(ESMF_GridComp)                   :: rof_comp

   type(ESMF_Array)                      :: x2da, d2xa, doma

   !----------------------------------------------------------
   
   !----------------------------------------------------------------------
   ! Determine cdata points
   !----------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
       gsMap=gsmap, dom=dom, infodata=infodata)

   call seq_comm_getcompstates(COMPID, rof_comp, import_state, export_state)
   call seq_infodata_GetData(infodata,rof_phase=phase)

   ! Copy infodata to state

   call esmfshr_infodata_infodata2state(infodata,export_state,ID=COMPID,rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase > 1) then
      call ESMF_StateGet(import_state, itemName="x2d", array=x2da, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call mct2esmf_copy(x2d, x2da, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

   ! call into ESMF init method
   call ESMF_GridCompInitialize(rof_comp, importState=import_state, exportState=export_state, &
        clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc,  endflag=ESMF_END_ABORT)

   ! copy export_state to infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase == 1) then
      call ESMF_StateGet(export_state, itemName="domain", array=doma, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(export_state, itemName="d2x", array=d2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(import_state, itemName="x2d", array=x2da, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(d2xa, COMPID, gsmap, mpicom, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(doma, dom, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(doma, dom%data, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(d2xa, d2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(d2xa, d2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
  
      call esmf2mct_init(x2da, x2d, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct_aVect_zero(x2d)
   else
      call ESMF_StateGet(export_state, itemName="d2x", array=d2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(d2xa, d2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

end subroutine rof_init_mct

!===============================================================================

subroutine rof_run_mct( EClock, cdata, x2d, d2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d
   type(mct_aVect)             ,intent(inout) :: d2x

   !----- local -----
   integer(IN)                      :: COMPID
   type(seq_infodata_type), pointer :: infodata
   type(ESMF_Array)                 :: d2xa,x2da 
   integer                          :: rc, urc
   type(ESMF_State)                 :: import_state, export_state
   type(ESMF_GridComp)              :: rof_comp
   !----------------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=COMPID, infodata=infodata)

   call seq_comm_getcompstates(COMPID, rof_comp, import_state, export_state)
  
   call esmfshr_infodata_infodata2state(infodata, export_state,ID=COMPID,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values to x2d
   call ESMF_StateGet(import_state, itemName="x2d", array=x2da, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call mct2esmf_copy(x2d, x2da, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! run rof component
   call ESMF_GridCompRun(rof_comp, importState=import_state, exportState=export_state, &
        clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc,  endflag=ESMF_END_ABORT)
   
   ! convert state back to infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values back to d2x
   call ESMF_StateGet(export_state, itemName="d2x", array=d2xa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call esmf2mct_copy(d2xa, d2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine rof_run_mct

!===============================================================================

subroutine rof_final_mct( EClock, cdata, x2d, d2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock  
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d
    type(mct_aVect)             ,intent(inout) :: d2x

    integer             :: rc, urc, COMPID
    type(ESMF_State)    :: import_state, export_state
    type(ESMF_GridComp) :: rof_comp

    !----------------------------------------------------------------------------
    ! Finalize routine 
    !----------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID)
    call seq_comm_getcompstates(COMPID, rof_comp, import_state, export_state)

    call ESMF_GridCompFinalize(rof_comp, importState=import_state, exportState=export_state, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc,  endflag=ESMF_END_ABORT)

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

