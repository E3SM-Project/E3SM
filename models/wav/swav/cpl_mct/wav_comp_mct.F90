module wav_comp_mct

! !USES:

  use mct_mod
  use esmf
  use seq_cdata_mod
  use seq_infodata_mod

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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_init_mct
!
! !DESCRIPTION:
!     stub wav model init
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine wav_init_mct( EClock, cdata, x2d, d2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename

!EOP
!-------------------------------------------------------------------------------

   call seq_infodata_PutData(cdata%infodata, wav_present=.false., &
        wav_prognostic=.false.)

end subroutine wav_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_run_mct
!
! !DESCRIPTION:
!     stub wav model run
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_run_mct( EClock, cdata, x2d, d2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        
   type(mct_aVect)             ,intent(inout) :: d2x        

!EOP
!-------------------------------------------------------------------------------

end subroutine wav_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_final_mct
!
! !DESCRIPTION:
!     stub wav model finalize 
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine wav_final_mct( EClock, cdata, x2d, d2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d
   type(mct_aVect)             ,intent(inout) :: d2x

!EOP
!-------------------------------------------------------------------------------

 end subroutine wav_final_mct

!===============================================================================

end module wav_comp_mct
