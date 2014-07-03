module wav_comp_esmf

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

  public :: wav_init_esmf
  public :: wav_run_esmf
  public :: wav_final_esmf
  public :: wav_register_esmf

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine wav_register_esmf(comp, rc)
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    print *, "In wav register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      wav_init_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      wav_run_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      wav_final_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)


end subroutine

!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_init_esmf
!
! !DESCRIPTION:
!     initialize dead wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_init_esmf(comp, import_state, export_state, EClock, rc)

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
    call ESMF_AttributeSet(export_state, name="wav_present", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="wav_prognostic", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "SWAV", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Sea Wav Stub Model", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Sea Wav", &
                           convention=convCIM, purpose=purpComp, rc=rc)

!    call ESMF_AttributeSet(comp, "Name", "someone", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "EmailAddress", &
!                           "someone@someplace", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
!                           convention=convCIM, purpose=purpComp, rc=rc)

end subroutine wav_init_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_run_esmf
!
! !DESCRIPTION:
!     run method for dead wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_run_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc


!EOP

   rc = ESMF_SUCCESS

end subroutine wav_run_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_final_esmf
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_final_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
 

   rc = ESMF_SUCCESS

end subroutine wav_final_esmf
!===============================================================================

end module wav_comp_esmf
