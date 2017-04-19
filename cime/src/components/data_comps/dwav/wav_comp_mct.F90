module wav_comp_mct

! !USES:

  use seq_cdata_mod
  use esmf
  use mct_mod

  use dwav_comp_mod

! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: wav_init_mct
  public :: wav_run_mct
  public :: wav_final_mct

!--------------------------------------------------------------------------
! Private data
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
!     initialize data wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_init_mct( EClock, cdata, x2w, w2x, NLFilename )

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2w, w2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    character(*), parameter :: subName = "(wav_init_mct) "
!-------------------------------------------------------------------------------


    if (present(NLFilename)) then
       call dwav_comp_init(EClock, cdata, x2w, w2x, NLFilename)
    else
       call dwav_comp_init(EClock, cdata, x2w, w2x)
    endif

end subroutine wav_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_run_mct
!
! !DESCRIPTION:
!     run method for data wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_run_mct( EClock, cdata,  x2w, w2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2w        ! driver -> data
   type(mct_aVect)             ,intent(inout) :: w2x        ! data   -> driver

!EOP

   character(*), parameter :: subName = "(wav_run_mct) "
!-------------------------------------------------------------------------------

   call dwav_comp_run(EClock, cdata, x2w, w2x)

end subroutine wav_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_final_mct
!
! !DESCRIPTION:
!     finalize method for data wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine wav_final_mct(EClock, cdata, x2w, w2x)

    implicit none

    !----- arguments -----

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2w        ! driver -> data
    type(mct_aVect)             ,intent(inout) :: w2x        ! data   -> driver

!EOP

    !--- formats ---
    character(*), parameter :: subName = "(wav_final_mct) "
!-------------------------------------------------------------------------------

    call dwav_comp_final()

end subroutine wav_final_mct
!===============================================================================
!===============================================================================


end module wav_comp_mct
