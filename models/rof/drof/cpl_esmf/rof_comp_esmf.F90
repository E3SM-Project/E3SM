module rof_comp_esmf

   use shr_kind_mod, only:  R8=>SHR_KIND_R8, IN=>SHR_KIND_IN, &
                            CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use shr_sys_mod   ! shared system calls

   use seq_cdata_mod
   use seq_infodata_mod

   use esmfshr_mod

   use drof_comp_mod
   use esmf
   use perf_mod
   use mct_mod

   implicit none

   public :: rof_init_esmf
   public :: rof_run_esmf
   public :: rof_final_esmf
   public :: rof_register_esmf

   private ! except

   type(seq_infodata_type)  :: infodata
   type(seq_cdata)     :: cdata
   type(mct_gsMap)     :: gsmap
   type(mct_gGrid)     :: ggrid
   type(mct_aVect)     :: x2d
   type(mct_aVect)     :: d2x

   !----- formats -----
   character(*),parameter :: subName =  "(rof_comp_esmf) "

   save ! save everything

!
! Author: Mariana Vertenstein
! This module is ESMF compliant rof data component 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

subroutine rof_register_esmf(comp, rc)

   implicit none

   type(ESMF_GridComp)  :: comp
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   print *, "In rof register routine"
   ! Register the callback routines.

   call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, rof_init_esmf, phase=1, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, rof_run_esmf, phase=1, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, rof_final_esmf, phase=1, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine

!===============================================================================

subroutine rof_init_esmf(comp, import_state, export_state, EClock, rc)
   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
   
   !----- local -----
   integer(IN)      :: MYID
   character(CL)    :: NLFilename
   type(ESMF_Array) :: Ex2d, Ed2x, Edom
   integer(IN)      :: phase
   
   character(*),parameter :: subName = "(rof_init_esmf) "
   character(ESMF_MAXSTR) :: convCIM, purpComp
   !----------------------------------------------------------
   
   rc = ESMF_SUCCESS

   NLFilename = 'unused'

   call esmfshr_infodata_state2infodata(export_state,infodata,ID=MYID,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call seq_infodata_GetData(infodata,rof_phase=phase)

   if (phase == 1) then
      call seq_cdata_init(cdata,MYID,ggrid,gsmap,infodata,'drof')
   else
      call ESMF_StateGet(import_state, itemName="x2d", array=Ex2d, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call esmf2mct_copy(Ex2d, x2d, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

   call drof_comp_init(EClock, cdata, x2d, d2x, NLFilename)

   call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   if (phase == 1) then
      Edom = mct2esmf_init(ggrid%data,gsmap,name='domain',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call mct2esmf_copy(ggrid%data,Edom,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      Ed2x = mct2esmf_init(d2x,gsmap,name='d2x',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call mct2esmf_copy(d2x,Ed2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      Ex2d = mct2esmf_init(x2d,gsmap,name='x2d',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateAdd(export_state,(/Edom/),rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateAdd(export_state,(/Ed2x/),rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call ESMF_StateAdd(import_state,(/Ex2d/),rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   else
      call ESMF_StateGet(export_state, itemName="d2x", array=Ed2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      call mct2esmf_copy(d2x,Ed2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   endif

    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "DROF", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Climatological River Runoff Data Model", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                 "The CESM data models perform the basic function of " // &
                 "reading external data, modifying that data, and then " // &
                 "sending it to the driver via standard CESM coupling " // &
                 "interfaces. The driver and other models have no " // &
                 "fundamental knowledge of whether another component " // &
                 "is fully active or just a data model.  In some cases, " // &
                 "data models are prognostic and also receive and use " // &
                 "some data sent by the driver to the data model.  But " // &
                 "in most cases, the data models are not running " // &
                 "prognostically and have no need to receive any data " // &
                 "from the driver.", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Runoff", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    !    call ESMF_AttributeSet(comp, "Name", "Sam Levis", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "EmailAddress", &
    !                           "slevis@ucar.edu", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    
   rc = ESMF_SUCCESS

end subroutine rof_init_esmf

!===============================================================================

subroutine rof_run_esmf(comp, import_state, export_state, EClock, rc)

   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

   !----- local -----
   integer(IN)      :: MYID
   type(ESMF_Array) :: Ex2d, Ed2x
   
   character(*),parameter :: subName = "(rof_run_esmf) "
   !----------------------------------------------------------
   
   rc = ESMF_SUCCESS

   ! Unpack import state

   call esmfshr_infodata_state2infodata(export_state,infodata,ID=MYID,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_StateGet(import_state, itemName="x2d", array=Ex2d, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call esmf2mct_copy(Ex2d, x2d, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! Run model

   call drof_comp_run(EClock, cdata, x2d, d2x)

   ! Pack export state

   call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_StateGet(export_state, itemName="d2x", array=Ed2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call mct2esmf_copy(d2x,Ed2x,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   rc = ESMF_SUCCESS

end subroutine rof_run_esmf

!===============================================================================

subroutine rof_final_esmf(comp, import_state, export_state, EClock, rc)

   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

   !----------------------------------------------------------------------------
   ! Finalize routine 
   !----------------------------------------------------------------------------
  
   rc = ESMF_SUCCESS

   call drof_comp_final()
  
end subroutine rof_final_esmf

!===============================================================================

end module rof_comp_esmf
