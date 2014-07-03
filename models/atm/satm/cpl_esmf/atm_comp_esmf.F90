module atm_comp_esmf

! !USES:

  use ESMF
  use esmfshr_mod
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_esmf
  public :: atm_run_esmf
  public :: atm_final_esmf
  public :: atm_register_esmf

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine atm_register_esmf(comp, rc)
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    print *, "In atm register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      atm_init_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      atm_run_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      atm_final_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)


end subroutine

!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_init_esmf
!
! !DESCRIPTION:
!     initialize dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_init_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

    ! Local variables
    character(ESMF_MAXSTR) :: convCIM, purpComp

!EOP

    rc = ESMF_SUCCESS

    ! Set flag to specify dead components
    call ESMF_AttributeSet(export_state, name="atm_present", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="atm_prognostic", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "SATM", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Atmosphere Stub Model", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Atmosphere", &
                           convention=convCIM, purpose=purpComp, rc=rc)

!    call ESMF_AttributeSet(comp, "IndividualName", "Cecile Hannay", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "IndividualEmailAddress", &
!                           "hannay@ucar.edu", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
!                           convention=convCIM, purpose=purpComp, rc=rc)

end subroutine atm_init_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_run_esmf
!
! !DESCRIPTION:
!     run method for dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_run_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc


!EOP

   rc = ESMF_SUCCESS

end subroutine atm_run_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_final_esmf
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_final_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
 

   rc = ESMF_SUCCESS

end subroutine atm_final_esmf
!===============================================================================

end module atm_comp_esmf
