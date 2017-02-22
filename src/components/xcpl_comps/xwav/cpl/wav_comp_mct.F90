module wav_comp_mct

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

  public :: wav_init_mct
  public :: wav_run_mct
  public :: wav_final_mct

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_init_mct
!
! !DESCRIPTION:
!     initialize dead wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine wav_init_mct( EClock, cdata, x2d, d2x, NLFilename )
     use seq_flds_mod     , only: flds_d2x => seq_flds_w2x_fields, &
         flds_x2d => seq_flds_x2w_fields

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP
    call dead_init_mct('wav', flds_d2x, flds_x2d, EClock, cdata, x2d, d2x, NLFilename )

end subroutine wav_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_run_mct
!
! !DESCRIPTION:
!     run method for dead wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_run_mct(EClock, cdata, x2d, d2x)

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver
!   type(eshr_timeMgr_clockType),intent(inout) :: SyncClock  ! Synchronization clock

!EOP

   call dead_run_mct('wav',EClock, cdata, x2d, d2x)


end subroutine wav_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_final_mct
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_final_mct(EClock, cdata, x2d, d2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

!EOP
    call dead_final_mct('wav', EClock, cdata, x2d, d2x)


end subroutine wav_final_mct
!===============================================================================

end module wav_comp_mct
