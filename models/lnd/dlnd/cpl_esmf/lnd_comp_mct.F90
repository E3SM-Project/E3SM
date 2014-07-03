!===============================================================================
! SVN: $Id: lnd_comp_mct.F90 16107 2009-05-18 17:39:58Z fei.liu@gmail.com $ 
! SVN: $URL: https://svn-ccsm-models.cgd.ucar.edu/dlnd7/branches/cpl7esmf_beta10/lnd_comp_mct.F90 $
!===============================================================================

module lnd_comp_mct

   use shr_kind_mod, only: CS=>shr_kind_CS, IN=>shr_kind_IN, R8=>shr_kind_R8

   use mct_mod
   use esmf
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_comm_mct     , only: seq_comm_getcompstates

   use esmfshr_mod
   use lnd_comp_esmf

   implicit none

   public :: lnd_init_mct
   public :: lnd_run_mct
   public :: lnd_final_mct
   public :: lnd_register

   private ! except

   save ! save everything

!
! Author: Fei Liu
! This module is a wrapper layer between ccsm driver and ESMF data lnd component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

subroutine lnd_register(lnd_petlist, ccsmComp, lnd_comp, import_state, export_state)

   implicit none

   integer, pointer                  :: lnd_petlist(:)
   type(ESMF_CplComp)                :: ccsmComp
   type(ESMF_GridComp),intent(out)   :: lnd_comp
   type(ESMF_State), intent(out)     :: import_state
   type(ESMF_State), intent(out)     :: export_state

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

  subroutine lnd_init_mct( EClock, cdata_l, x2l, l2x, cdata_s, x2s, s2x, NLFilename )

   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_Clock),intent(inout)              :: EClock
   type(seq_cdata), intent(inout)              :: cdata_l
   type(mct_aVect), intent(inout)              :: x2l
   type(mct_aVect), intent(inout)              :: l2x   
   type(seq_cdata), intent(inout)              :: cdata_s
   type(mct_aVect), intent(inout)              :: x2s
   type(mct_aVect), intent(inout)              :: s2x   
   character(len=*), optional,   intent(in)    :: NLFilename ! Namelist filename
   
   !----- local -----
   integer(IN)                           :: MYID
   integer(IN)                           :: mpicom
   type(mct_gsMap)             , pointer :: gsmap_l
   type(mct_gGrid)             , pointer :: dom_l
   type(mct_gsMap)             , pointer :: gsmap_s
   type(mct_gGrid)             , pointer :: dom_s
   type(seq_infodata_type), pointer      :: infodata
   integer                               :: rc, urc
   integer(IN)                           :: phase
   type(ESMF_State)                      :: import_state, export_state
   type(ESMF_GridComp)                   :: lnd_comp

   type(ESMF_Array)                      :: x2la, l2xa, domla
   type(ESMF_Array)                      ::       r2xa, domra
   type(ESMF_Array)                      :: x2sa, s2xa, domsa

   character(*),parameter :: subName = "(lnd_init_mct) "
   !----------------------------------------------------------
   
   !----------------------------------------------------------------------
   ! Determine cdata points
   !----------------------------------------------------------------------

   call seq_cdata_setptrs(cdata_l, ID=MYID, mpicom=mpicom, &
       gsMap=gsmap_l, dom=dom_l, infodata=infodata)
   call seq_cdata_setptrs(cdata_s, gsMap=gsmap_s, dom=dom_s)

   call seq_comm_getcompstates(MYID, lnd_comp, import_state, export_state)
   call seq_infodata_GetData(infodata,lnd_phase=phase)

   ! Copy infodata to state

   call esmfshr_infodata_infodata2state(infodata,export_state,ID=MYID,rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase > 1) then
      call ESMF_StateGet(import_state, itemName="x2l", array=x2la, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct2esmf_copy(x2l, x2la, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call ESMF_StateGet(import_state, itemName="x2s", array=x2sa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct2esmf_copy(x2s, x2sa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

   ! call into ESMF init method
   call ESMF_GridCompInitialize(lnd_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy export_state to infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase == 1) then
      call ESMF_StateGet(export_state, itemName="domain_l", array=domla, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(export_state, itemName="l2x", array=l2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(import_state, itemName="x2l", array=x2la, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(export_state, itemName="domain_s", array=domsa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(export_state, itemName="s2x", array=s2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(import_state, itemName="x2s", array=x2sa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(l2xa, MYID, gsmap_l, mpicom, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(s2xa, MYID, gsmap_s, mpicom, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(domla, dom_l, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(domla, dom_l%data, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(domsa, dom_s, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(domsa, dom_s%data, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(l2xa, l2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(l2xa, l2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(s2xa, s2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(s2xa, s2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(x2la, x2l, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct_aVect_zero(x2l)

      call esmf2mct_init(x2sa, x2s, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct_aVect_zero(x2s)
   else
      call ESMF_StateGet(export_state, itemName="l2x", array=l2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(l2xa, l2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(export_state, itemName="s2x", array=s2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(s2xa, s2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

end subroutine lnd_init_mct

!===============================================================================

subroutine lnd_run_mct( EClock, cdata_l, x2l, l2x, cdata_s, x2s, s2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata_l
   type(mct_aVect)             ,intent(inout) :: x2l
   type(mct_aVect)             ,intent(inout) :: l2x
   type(seq_cdata)             ,intent(inout) :: cdata_s
   type(mct_aVect)             ,intent(inout) :: x2s
   type(mct_aVect)             ,intent(inout) :: s2x

   !----- local -----
   integer(IN)                      :: MYID
   type(seq_infodata_type), pointer :: infodata
   type(ESMF_Array)                 :: l2xa,x2la
   type(ESMF_Array)                 :: s2xa,x2sa
   integer                          :: rc, urc
   type(ESMF_State)                 :: import_state, export_state
   type(ESMF_GridComp)              :: lnd_comp
   !----------------------------------------------------------------------------

   call seq_cdata_setptrs(cdata_l, ID=MYID, infodata=infodata)

   call seq_comm_getcompstates(MYID, lnd_comp, import_state, export_state)
  
   call esmfshr_infodata_infodata2state(infodata, export_state,ID=MYID,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values to x2l
   call ESMF_StateGet(import_state, itemName="x2l", array=x2la, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call mct2esmf_copy(x2l, x2la, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_StateGet(import_state, itemName="x2s", array=x2sa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call mct2esmf_copy(x2s, x2sa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_GridCompRun(lnd_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc,  endflag=ESMF_END_ABORT)
   
   ! convert state back to infodata, the new nextsw_cday is updated in infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values back to l2x
   call ESMF_StateGet(export_state, itemName="l2x", array=l2xa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call esmf2mct_copy(l2xa, l2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_StateGet(export_state, itemName="s2x", array=s2xa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call esmf2mct_copy(s2xa, s2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine lnd_run_mct

!===============================================================================

subroutine lnd_final_mct( EClock, cdata, x2d, d2x, cdata_s, x2s, s2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock    
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d       
    type(mct_aVect)             ,intent(inout) :: d2x       
    type(seq_cdata)             ,intent(inout) :: cdata_s
    type(mct_aVect)             ,intent(inout) :: x2s
    type(mct_aVect)             ,intent(inout) :: s2x

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
    if (rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc,  endflag=ESMF_END_ABORT)

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

