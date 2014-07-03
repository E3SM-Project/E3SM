module rof_comp_mct

   use shr_kind_mod, only: CS=>shr_kind_CS, IN=>shr_kind_IN, R8=>shr_kind_R8
   use perf_mod    , only : t_startf, t_stopf, t_barrierf

   use mct_mod
   use esmf
   use seq_comm_mct, only: seq_comm_getcompstates
   use seq_cdata_mod
   use seq_infodata_mod
   use shr_sys_mod

   use esmfshr_mod
   use rof_comp_esmf
   use RtmVar, only : iulog

   implicit none

   public :: rof_init_mct
   public :: rof_run_mct
   public :: rof_final_mct
   public :: rof_register

   private ! except

   save ! save everything
!
! Author: Mariana Vertenstein
! This module is a wrapper layer between cesm driver and ESMF rof  component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

  subroutine rof_register(rof_petlist, ccsmComp, rof_comp, import_state, export_state)

    implicit none

    integer, pointer                  :: rof_petlist(:)
    type(ESMF_CplComp) ,intent(inout) :: ccsmComp
    type(ESMF_GridComp),intent(out)   :: rof_comp
    type(ESMF_State)   ,intent(out)   :: import_state, export_state

    integer :: rc

    rof_comp = ESMF_GridCompCreate(name="rof_comp", petList=rof_petlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create rof comp')

    call ESMF_GridCompSetServices(rof_comp, rof_register_esmf, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register rof comp')

    import_state = ESMF_StateCreate(name="rof import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import rof state')

    export_state = ESMF_StateCreate(name="rof export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export rof state')

    call ESMF_AttributeLink(ccsmComp, rof_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

  end subroutine rof_register

!===============================================================================

  subroutine rof_init_mct( EClock, cdata, x2r, r2x, NLFilename )

   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_Clock),intent(inout)              :: EClock
   type(seq_cdata), intent(inout)              :: cdata
   type(mct_aVect), intent(inout)              :: x2r   
   type(mct_aVect), intent(inout)              :: r2x   
   character(len=*), optional,   intent(in)    :: NLFilename ! Namelist filename
   
   !----- local -----
   integer(IN)                           :: ROFID
   integer(IN)                           :: mpicom
   type(mct_gsMap)             , pointer :: gsmap
   type(mct_gGrid)             , pointer :: dom
   type(seq_infodata_type), pointer      :: infodata
   integer(IN)                           :: gsize
   integer                               :: rc, urc
   integer(IN)                           :: phase
   logical                               :: rof_present

   type(ESMF_Array)                      :: x2ra, r2xa, doma
   type(ESMF_State)                      :: import_state, export_state
   type(ESMF_GridComp)                   :: rof_comp

   character(*),parameter :: subName = "(rof_init_mct) "
   !----------------------------------------------------------
   
   ! Determine cdata points

   call seq_cdata_setptrs(cdata, ID=ROFID, mpicom=mpicom, &
        gsMap=gsmap, dom=dom, infodata=infodata)
   call seq_cdata_setptrs(cdata, gsMap=gsmap, dom=dom)
   call seq_comm_getcompstates(ROFID, rof_comp, import_state, export_state)
   call seq_infodata_GetData(infodata, rof_phase=phase)

   ! Copy infodata to state

   call esmfshr_infodata_infodata2state(infodata,export_state,ID=ROFID,rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase > 1) then
      call ESMF_StateGet(import_state, itemName="x2r", array=x2ra, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call mct2esmf_copy(x2r, x2ra, rc=rc)
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

   call seq_infodata_GetData(infodata, rof_present=rof_present)
   if (.not. rof_present) return

   if (phase == 1) then
      call ESMF_StateGet(export_state, itemName="domain", array=doma, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(export_state, itemName="r2x", array=r2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(import_state, itemName="x2r", array=x2ra, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_AttributeGet(export_state, name="gsize", value=gsize, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(r2xa, ROFID, gsmap, mpicom, gsize=gsize, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(doma, dom, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call esmf2mct_copy(doma, dom%data, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(r2xa, r2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call esmf2mct_copy(r2xa, r2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
  
      call esmf2mct_init(x2ra, x2r, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct_aVect_zero(x2r)
   else
      call ESMF_StateGet(export_state, itemName="r2x", array=r2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call esmf2mct_copy(r2xa, r2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

end subroutine rof_init_mct

!===============================================================================

subroutine rof_run_mct( EClock, cdata, x2r, r2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)  ,intent(inout) :: EClock
   type(seq_cdata)   ,intent(inout) :: cdata
   type(mct_aVect)   ,intent(inout) :: x2r
   type(mct_aVect)   ,intent(inout) :: r2x

   !----- local -----
   type(seq_infodata_type), pointer :: infodata
   type(ESMF_Array)                 :: r2xa,x2ra,doma  
   type(mct_gGrid), pointer         :: dom
   integer(IN)                      :: ROFID
   integer(IN)                      :: rc, urc
   type(ESMF_State)                 :: import_state, export_state
   type(ESMF_GridComp)              :: rof_comp
   logical                          :: rof_present
   !----------------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=ROFID, infodata=infodata)
   call seq_infodata_GetData(infodata, rof_present=rof_present)
   if (.not. rof_present) return

   call seq_comm_getcompstates(ROFID, rof_comp, import_state, export_state)

   call t_startf('rtm_run1')
   call esmfshr_infodata_infodata2state(infodata, export_state, ID=ROFID, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values to x2r
   call ESMF_StateGet(import_state, itemName="x2r", array=x2ra, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call mct2esmf_copy(x2r, x2ra, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call t_stopf('rtm_run1')

   call t_startf('rtm_gcrun')
   call ESMF_GridCompRun(rof_comp, importState=import_state, exportState=export_state, &
        clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc,  endflag=ESMF_END_ABORT)
   call t_stopf('rtm_gcrun')

   call t_startf('rtm_run2')
   ! convert state back to infodata, the new nextsw_cday is updated in infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values back to r2x
   call ESMF_StateGet(export_state, itemName="r2x", array=r2xa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call esmf2mct_copy(r2xa, r2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call t_stopf('rtm_run2')

end subroutine rof_run_mct

!===============================================================================

subroutine rof_final_mct( EClock, cdata, x2r, r2x)

    implicit none

    !----- arguments -----

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2r
    type(mct_aVect)             ,intent(inout) :: r2x

    !----- local -----
    integer                                    :: rc, urc, COMPID
    type(seq_infodata_type), pointer           :: infodata
    type(ESMF_State)                           :: import_state, export_state
    type(ESMF_GridComp)                        :: rof_comp
    logical                                    :: rof_present

    !----------------------------------------------------------------------------
    ! Finalize routine 
    !----------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID, infodata=infodata)
    call seq_infodata_GetData(infodata, rof_present=rof_present)
    if (.not. rof_present) return

    call seq_comm_getcompstates(COMPID, rof_comp, import_state, export_state)

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

