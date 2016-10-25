module lnd_comp_esmf

#ifdef ESMF_INTERFACE 
  use shr_kind_mod, only:  R8=>SHR_KIND_R8, IN=>SHR_KIND_IN, &
       CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_sys_mod   ! shared system calls

  use seq_cdata_mod
  use seq_infodata_mod

  use esmf
  use esmfshr_mod

  use dlnd_comp_mod
  use perf_mod
  use mct_mod

  implicit none

  public :: lnd_init_esmf
  public :: lnd_run_esmf
  public :: lnd_final_esmf
  public :: lnd_register_esmf

  private ! except

  type(seq_infodata_type)  :: infodata
  type(seq_cdata)     :: cdata
  type(mct_gsMap)     :: gsmap
  type(mct_gGrid)     :: ggrid
  type(mct_aVect)     :: x2l
  type(mct_aVect)     :: l2x

  !----- formats -----
  character(*),parameter :: subName =  "(lnd_comp_esmf) "

  save ! save everything

  !
  ! Author: Fei Liu
  ! This module is ESMF compliant lnd data component 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================

  subroutine lnd_register_esmf(comp, rc)

    implicit none

    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    print *, "In lnd register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, lnd_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, lnd_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, lnd_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine lnd_register_esmf

  !===============================================================================

  subroutine lnd_init_esmf(comp, import_state, export_state, EClock, rc)
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
    type(ESMF_Array) :: Ex2l, El2x, Edoml

    character(*),parameter :: subName = "(lnd_init_esmf) "
    character(ESMF_MAXSTR) :: convCIM, purpComp
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    NLFilename = 'unused'

    call esmfshr_infodata_state2infodata(export_state,infodata,ID=MYID)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call seq_cdata_init(cdata, MYID, ggrid, gsmap, infodata,'dlnd')

    call dlnd_comp_init(EClock, cdata, x2l, l2x, NLFilename)

    call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    Edoml = mct2esmf_init(ggrid%data,gsmap,name='domain',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call mct2esmf_copy(ggrid%data,Edoml,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    El2x = mct2esmf_init(l2x,gsmap,name='d2x',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call mct2esmf_copy(l2x,El2x,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    Ex2l = mct2esmf_init(x2l,gsmap,name='x2d',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state,(/Edoml/),rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state,(/El2x/),rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(import_state,(/Ex2l/),rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)


#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
         convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "DLND", &
         convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
         "Climatological Land Data Model", &
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
    call ESMF_AttributeSet(comp, "ModelType", "Land", &
         convention=convCIM, purpose=purpComp, rc=rc)

    !    call ESMF_AttributeSet(comp, "Name", "Sam Levis", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "EmailAddress", &
    !                           "slevis@ucar.edu", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
#endif

    rc = ESMF_SUCCESS

  end subroutine lnd_init_esmf

  !===============================================================================

  subroutine lnd_run_esmf(comp, import_state, export_state, EClock, rc)

    implicit none

    !----- arguments -----
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

    !----- local -----
    integer(IN)      :: MYID
    type(ESMF_Array) :: Ex2l, El2x

    character(*),parameter :: subName = "(lnd_run_esmf) "
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Unpack import state

    call esmfshr_infodata_state2infodata(export_state,infodata,ID=MYID)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(import_state, itemName="x2d", array=Ex2l, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_copy(Ex2l, x2l, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Run model

    call dlnd_comp_run(EClock, cdata, x2l, l2x)

    ! Pack export state

    call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName="d2x", array=El2x, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call mct2esmf_copy(l2x,El2x,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  end subroutine lnd_run_esmf

  !===============================================================================

  subroutine lnd_final_esmf(comp, import_state, export_state, EClock, rc)

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

    call dlnd_comp_final()

  end subroutine lnd_final_esmf

  !===============================================================================
#endif

end module lnd_comp_esmf
