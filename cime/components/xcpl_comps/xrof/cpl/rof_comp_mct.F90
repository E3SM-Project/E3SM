module rof_comp_mct

! !USES:

  use dead_mct_mod, only : dead_init_mct, dead_run_mct, dead_final_mct

  use esmf, only : esmf_clock
  use seq_cdata_mod, only : seq_cdata
  use mct_mod, only : mct_avect

!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: rof_init_mct
  public :: rof_run_mct
  public :: rof_final_mct

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_init_mct
!
! !DESCRIPTION:
!     initialize dead rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine rof_init_mct( EClock, cdata, x2d, d2x, NLFilename )
    use seq_flds_mod     , only: flds_d2x => seq_flds_r2x_fields, &
         flds_x2d => seq_flds_x2r_fields

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP
    call dead_init_mct('rof', flds_d2x, flds_x2d, EClock, cdata, x2d, d2x, NLFilename )

end subroutine rof_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_run_mct
!
! !DESCRIPTION:
!     run method for dead rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine rof_run_mct(EClock, cdata, x2d, d2x)

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver
!   type(eshr_timeMgr_clockType),intent(inout) :: SyncClock  ! Synchronization clock

!EOP

   call dead_run_mct('rof',EClock, cdata, x2d, d2x)


end subroutine rof_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_final_mct
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine rof_final_mct(EClock, cdata, x2d, d2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

!EOP
    call dead_final_mct('rof', EClock, cdata, x2d, d2x)


end subroutine rof_final_mct
!===============================================================================

end module rof_comp_mct
