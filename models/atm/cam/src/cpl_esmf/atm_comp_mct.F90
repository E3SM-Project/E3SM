!===============================================================================
! SVN: $Id: atm_comp_mct.F90 16107 2009-05-18 17:39:58Z fei.liu@gmail.com $ 
! SVN: $URL: https://svn-ccsm-models.cgd.ucar.edu/datm7/branches/cpl7esmf_beta10/atm_comp_mct.F90 $
!===============================================================================

module atm_comp_mct

   use shr_kind_mod, only: CS=>shr_kind_CS, IN=>shr_kind_IN, R8=>shr_kind_R8

   use mct_mod
   use esmf
   use seq_comm_mct     , only: seq_comm_getcompstates
   use seq_cdata_mod
   use seq_infodata_mod

   use esmfshr_mod
   use atm_comp_esmf

   implicit none

   public :: atm_init_mct
   public :: atm_run_mct
   public :: atm_final_mct
   public :: atm_register

   private ! except

   type(ESMF_GridComp)     :: atm_comp
   type(ESMF_State)        :: import_state, export_state

   save ! save everything

!
! Author: Fei Liu
! This module is a wrapper layer between ccsm driver and ESMF data atm component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

subroutine atm_register(atm_petlist, ccsmComp, atm_comp, import_state, export_state)

    implicit none

    ! Arguments
    integer, pointer                  :: atm_petlist(:)
    type(ESMF_CplComp) ,intent(inout) :: ccsmComp
    type(ESMF_GridComp),intent(out)   :: atm_comp
    type(ESMF_State)   ,intent(out)   :: import_state, export_state

    ! Local variables
    integer                           :: rc

    atm_comp = ESMF_GridCompCreate(name="atm_comp", petList=atm_petlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create atm comp')
    call ESMF_GridCompSetServices(atm_comp, atm_register_esmf, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register atm comp')
    import_state = ESMF_StateCreate(name="atm import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import atm state')
    export_state = ESMF_StateCreate(name="atm export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export atm state')

    ! Link attribute
    call ESMF_AttributeLink(ccsmComp, atm_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

end subroutine

!===============================================================================

  subroutine atm_init_mct( EClock, cdata, x2d, d2x, NLFilename )

   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_Clock),intent(inout)              :: EClock
   type(seq_cdata), intent(inout)              :: cdata
   type(mct_aVect), intent(inout)              :: x2d
   type(mct_aVect), intent(inout)              :: d2x   
   character(len=*), optional,   intent(in)    :: NLFilename ! Namelist filename
   
   !----- local -----
   integer(IN)                           :: ATMID
   integer(IN)                           :: mpicom
   type(mct_gsMap)             , pointer :: gsMap
   type(mct_gGrid)             , pointer :: dom
   type(seq_infodata_type), pointer      :: infodata
   integer                               :: rc, urc
   integer(IN)                           :: phase

   type(ESMF_Array)                      :: x2da, d2xa, doma
   type(ESMF_State)                      :: import_state, export_state
   type(ESMF_GridComp)                   :: atm_comp

   character(*),parameter :: subName = "(atm_init_mct) "
   !----------------------------------------------------------
   
   !----------------------------------------------------------------------
   ! Determine cdata points
   !----------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=ATMID, mpicom=mpicom, &
       gsMap=gsMap, dom=dom, infodata=infodata)
   call seq_comm_getcompstates(ATMID, atm_comp, import_state, export_state)

   call seq_infodata_GetData(infodata,atm_phase=phase)

   ! Copy infodata to state

   call esmfshr_infodata_infodata2state(infodata, export_state, rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase > 1) then
      call ESMF_StateGet(import_state, itemName="x2d", array=x2da, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct2esmf_copy(x2d, x2da, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

   ! call into ESMF init method
   call ESMF_AttributeSet(export_state, name="ID", value=ATMID, rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_GridCompInitialize(atm_comp, importState=import_state, exportState=export_state, &
        clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

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

      ! Initialize MCT gsMap 
      call esmf2mct_init(d2xa, ATMID, gsMap, mpicom, rc=rc)
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

end subroutine atm_init_mct

!===============================================================================

subroutine atm_run_mct( EClock, cdata, x2d, d2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d
   type(mct_aVect)             ,intent(inout) :: d2x

   !----- local -----
   type(seq_infodata_type), pointer :: infodata
   type(ESMF_Array)                 :: d2xa,x2da
   integer                          :: rc, urc, COMPID
   type(ESMF_State)                 :: import_state, export_state
   type(ESMF_GridComp)              :: atm_comp
   integer                          :: ATMID
   !----------------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=COMPID, infodata=infodata)
   call seq_comm_getcompstates(COMPID, atm_comp, import_state, export_state)

   call esmfshr_infodata_infodata2state(infodata, export_state, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values to x2d
   call ESMF_StateGet(import_state, itemName="x2d", array=x2da, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call mct2esmf_copy(x2d, x2da, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(export_state, name="ID", value=ATMID, rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_GridCompRun(atm_comp, importState=import_state, exportState=export_state, &
        clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   
   ! convert state back to infodata, the new nextsw_cday is updated in infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values back to d2x
   call ESMF_StateGet(export_state, itemName="d2x", array=d2xa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call esmf2mct_copy(d2xa, d2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine atm_run_mct

!===============================================================================

subroutine atm_final_mct(EClock, cdata, x2d, d2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

   !----- local -----
   integer                                    :: rc, urc
   integer                                    :: COMPID
   type(ESMF_State)                           :: import_state, export_state
   type(ESMF_GridComp)                        :: atm_comp
   !----------------------------------------------------------------------------
   ! Finalize routine 
   !----------------------------------------------------------------------------
   call seq_cdata_setptrs(cdata, ID=COMPID)
   call seq_comm_getcompstates(COMPID, atm_comp, import_state, export_state)
    
   call ESMF_GridCompFinalize(atm_comp, importState=import_state, exportState=export_state, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! destroy component and states
   call ESMF_StateDestroy(import_state, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_StateDestroy(export_state, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_GridCompDestroy(atm_comp, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine atm_final_mct

!===============================================================================

end module atm_comp_mct

