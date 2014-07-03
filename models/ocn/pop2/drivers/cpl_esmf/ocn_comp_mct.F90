!===============================================================================
! SVN: $Id: ocn_comp_mct.F90 16107 2009-05-18 17:39:58Z fei.liu@gmail.com $ 
! SVN: $URL: https://svn-ccsm-models.cgd.ucar.edu/docn7/branches/cpl7esmf_beta10/ocn_comp_mct.F90 $
!===============================================================================

module ocn_comp_mct

   use shr_kind_mod, only: CS=>shr_kind_CS, IN=>shr_kind_IN, R8=>shr_kind_R8

   use mct_mod
   use esmf
   use seq_comm_mct     , only: seq_comm_getcompstates
   use seq_cdata_mod
   use seq_infodata_mod
   use esmfshr_mod
   use ocn_comp_esmf

   use POP_MCT_vars_mod

   implicit none

   public :: ocn_init_mct
   public :: ocn_run_mct
   public :: ocn_final_mct
   public :: ocn_register

   private ! except

   save ! save everything

!
! Author: Fei Liu
! This module is a wrapper layer between ccsm driver and ESMF data ocn component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

subroutine ocn_register(ocn_petlist, ccsmComp, ocn_comp, import_state, export_state)

    implicit none

    integer, pointer                  :: ocn_petlist(:)
    type(ESMF_CplComp) ,intent(inout) :: ccsmComp
    type(ESMF_GridComp),intent(out)   :: ocn_comp
    type(ESMF_State)   ,intent(out)   :: import_state, export_state

    integer            :: rc

    ocn_comp = ESMF_GridCompCreate(name="ocn_comp", petList=ocn_petlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create ocn comp')
    call ESMF_GridCompSetServices(ocn_comp, ocn_register_esmf, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to register ocn comp')
    import_state = ESMF_StateCreate(name="ocn import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create import ocn state')
    export_state = ESMF_StateCreate(name="ocn export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export ocn state')

    call ESMF_AttributeLink(ccsmComp, ocn_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attribute')

end subroutine

!===============================================================================

  subroutine ocn_init_mct( EClock, cdata, x2d, d2x, NLFilename )

   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_Clock),intent(inout)              :: EClock
   type(seq_cdata), intent(inout)              :: cdata
   type(mct_aVect), intent(inout)              :: x2d
   type(mct_aVect), intent(inout)              :: d2x   
   character(len=*), optional,   intent(in)    :: NLFilename ! Namelist filename
   
   !----- local -----
   integer(IN)                           :: OCNID
   integer(IN)                           :: mpicom
   type(mct_gsMap)             , pointer :: gsMap
   type(mct_gGrid)             , pointer :: dom
   type(seq_infodata_type), pointer      :: infodata
   integer(IN)                           :: gsize
   integer                               :: rc, urc
   integer(IN)                           :: phase

   type(ESMF_Array)                      :: x2da, d2xa, doma
   type(ESMF_State)                      :: import_state, export_state
   type(ESMF_GridComp)                   :: ocn_comp

   character(*),parameter :: subName = "(ocn_init_mct) "
   !----------------------------------------------------------
   
   !----------------------------------------------------------------------
   ! Determine cdata points
   !----------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=OCNID, mpicom=mpicom, &
       gsMap=gsMap, dom=dom, infodata=infodata)
   call seq_comm_getcompstates(OCNID, ocn_comp, import_state, export_state)

   POP_MCT_OCNID   =  OCNID
   POP_MCT_gsMap_o => gsMap
   POP_MCT_dom_o   => dom

   call seq_infodata_GetData(infodata,ocn_phase=phase)

   ! Copy infodata to state

   call esmfshr_infodata_infodata2state(infodata,export_state,id=OCNID,rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase > 1) then
      call ESMF_StateGet(import_state, itemName="x2d", array=x2da, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct2esmf_copy(x2d, x2da, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

   ! call into ESMF init method
   call ESMF_GridCompInitialize(ocn_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
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

      call ESMF_AttributeGet(export_state, name="gsize", value=gsize, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      ! Initialize MCT gsMap 
      call esmf2mct_init(d2xa, OCNID, gsMap, mpicom, gsize=gsize, rc=rc)
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

end subroutine ocn_init_mct

!===============================================================================

subroutine ocn_run_mct( EClock, cdata, x2d, d2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d
   type(mct_aVect)             ,intent(inout) :: d2x

   !----- local -----
   type(seq_infodata_type), pointer :: infodata
   type(ESMF_Array)                 :: d2xa,x2da
   integer                          :: rc, urc
   type(ESMF_State)                 :: import_state, export_state
   type(ESMF_GridComp)              :: ocn_comp
   integer                          :: OCNID
   !----------------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, id=OCNID, infodata=infodata)
   call seq_comm_getcompstates(OCNID, ocn_comp, import_state, export_state)

   call esmfshr_infodata_infodata2state(infodata, export_state, id=OCNID, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values to x2d
   call ESMF_StateGet(import_state, itemName="x2d", array=x2da, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call mct2esmf_copy(x2d, x2da, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_GridCompRun(ocn_comp, importState=import_state, exportState=export_state, clock=EClock, userRc=urc, rc=rc)
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

end subroutine ocn_run_mct

!===============================================================================

subroutine ocn_final_mct(EClock, cdata, x2d, d2x)

    implicit none

    !----- arguments -----
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    integer                               :: rc, urc, COMPID
    type(ESMF_State)                      :: import_state, export_state
    type(ESMF_GridComp)                   :: ocn_comp

    !----------------------------------------------------------------------------
    ! Finalize routine 
    !----------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID)
    call seq_comm_getcompstates(COMPID, ocn_comp, import_state, export_state)

    call ESMF_GridCompFinalize(ocn_comp, importState=import_state, exportState=export_state, userRc=urc, rc=rc)
    if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! destroy component and states
    call ESMF_StateDestroy(import_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateDestroy(export_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompDestroy(ocn_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine ocn_final_mct

!===============================================================================

end module ocn_comp_mct

