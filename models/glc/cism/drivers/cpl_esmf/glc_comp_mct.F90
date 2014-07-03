!===============================================================================
! SVN: $Id: ice_comp_mct.F90 16107 2009-05-18 17:39:58Z fei.liu@gmail.com $ 
! SVN: $URL: https://svn-ccsm-models.cgd.ucar.edu/dice7/branches/cpl7esmf_beta10/ice_comp_mct.F90 $
!===============================================================================

module glc_comp_mct

   use shr_kind_mod, only: CS=>shr_kind_CS, IN=>shr_kind_IN, R8=>shr_kind_R8

   use mct_mod
   use esmf
   use seq_comm_mct     , only: seq_comm_getcompstates
   use seq_cdata_mod
   use seq_infodata_mod

   use esmfshr_mod
   use glc_comp_esmf

   implicit none

   public :: glc_init_mct
   public :: glc_run_mct
   public :: glc_final_mct
   public :: glc_register

   private ! except

   type(ESMF_GridComp)     :: glc_comp
   type(ESMF_State)        :: import_state, export_state

   save ! save everything

!
! Author: Fei Liu
! This module is a wrapper layer between ccsm driver and ESMF data glc component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

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
    if(rc /= ESMF_SUCCESS) call shr_sys_abort('failed to create export glc state')

end subroutine

!===============================================================================

  subroutine glc_init_mct( EClock, cdata, x2g, g2x, NLFilename )

   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_Clock),intent(inout)              :: EClock
   type(seq_cdata), intent(inout)              :: cdata
   type(mct_aVect), intent(inout)              :: x2g
   type(mct_aVect), intent(inout)              :: g2x   
   character(len=*), optional,   intent(in)    :: NLFilename ! Namelist filename
   
   !----- local -----
   integer(IN)                           :: GLCID
   integer(IN)                           :: mpicom
   type(mct_gsMap)             , pointer :: gsMap
   type(mct_gGrid)             , pointer :: dom
   type(seq_infodata_type), pointer      :: infodata
   integer(IN)                           :: gsize
   integer                               :: rc, urc
   integer(IN)                           :: phase

   type(ESMF_Array)                      :: x2ga, g2xa, doma
   type(ESMF_State)                      :: import_state, export_state
   type(ESMF_GridComp)                   :: glc_comp

   character(*),parameter :: subName = "(glc_init_mct) "
   !----------------------------------------------------------
   
   !----------------------------------------------------------------------
   ! Determine cdata points
   !----------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=GLCID, mpicom=mpicom, &
       gsMap=gsMap, dom=dom, infodata=infodata)
   call seq_comm_getcompstates(GLCID, glc_comp, import_state, export_state)

   call seq_infodata_GetData(infodata,glc_phase=phase)

   ! Copy infodata to state

   call esmfshr_infodata_infodata2state(infodata,export_state,ID=GLCID,rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase > 1) then
      call ESMF_StateGet(import_state, itemName="x2g", array=x2ga, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct2esmf_copy(x2g, x2ga, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

   ! call into ESMF init method
   call ESMF_GridCompInitialize(glc_comp, importState=import_state, exportState=export_state, &
        clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy export_state to infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase == 1) then
      call ESMF_StateGet(export_state, itemName="domain", array=doma, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(export_state, itemName="g2x", array=g2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateGet(import_state, itemName="x2g", array=x2ga, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_AttributeGet(export_state, name="gsize", value=gsize, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      ! Initialize MCT gsMap 
      call esmf2mct_init(g2xa, GLCID, gsMap, mpicom, gsize=gsize, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(doma, dom, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(doma, dom%data, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_init(g2xa, g2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(g2xa, g2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
  
      call esmf2mct_init(x2ga, x2g, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct_aVect_zero(x2g)
   else
      call ESMF_StateGet(export_state, itemName="g2x", array=g2xa, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(g2xa, g2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

end subroutine glc_init_mct

!===============================================================================

subroutine glc_run_mct( EClock, cdata, x2g, g2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2g
   type(mct_aVect)             ,intent(inout) :: g2x

   !----- local -----
   type(seq_infodata_type), pointer :: infodata
   type(ESMF_Array)                 :: g2xa,x2ga
   integer                          :: rc, urc, COMPID
   type(ESMF_State)                 :: import_state, export_state
   type(ESMF_GridComp)              :: glc_comp
   integer                          :: GLCID
   !----------------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=COMPID, infodata=infodata)
   call seq_comm_getcompstates(COMPID, glc_comp, import_state, export_state)

   call esmfshr_infodata_infodata2state(infodata, export_state, ID=COMPID, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values to x2g
   call ESMF_StateGet(import_state, itemName="x2g", array=x2ga, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call mct2esmf_copy(x2g, x2ga, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_GridCompRun(glc_comp, importState=import_state, exportState=export_state, &
        clock=EClock, userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   
   ! convert state back to infodata, the new nextsw_cday is updated in infodata
   call esmfshr_infodata_state2infodata(export_state, infodata, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! copy values back to g2x
   call ESMF_StateGet(export_state, itemName="g2x", array=g2xa, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call esmf2mct_copy(g2xa, g2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine glc_run_mct

!===============================================================================

subroutine glc_final_mct( EClock, cdata, x2g, g2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2g 
   type(mct_aVect)             ,intent(inout) :: g2x 

   !----- local -----
   integer                          :: rc, urc, COMPID
   type(ESMF_State)                 :: import_state, export_state
   type(ESMF_GridComp)              :: glc_comp
   !----------------------------------------------------------------------------

   call seq_cdata_setptrs(cdata, ID=COMPID)
   call seq_comm_getcompstates(COMPID, glc_comp, import_state, export_state)

   !----------------------------------------------------------------------------
   ! Finalize routine 
   !----------------------------------------------------------------------------
    
   call ESMF_GridCompFinalize(glc_comp, importState=import_state, exportState=export_state, &
        userRc=urc, rc=rc)
   if(urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

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

