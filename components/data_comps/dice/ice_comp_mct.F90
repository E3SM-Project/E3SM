module ice_comp_mct

! !USES:

  use seq_cdata_mod
  use esmf
  use mct_mod

  use dice_comp_mod

! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: ice_init_mct
  public :: ice_run_mct
  public :: ice_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_init_mct
!
! !DESCRIPTION:
!     initialize data ice model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ice_init_mct( EClock, cdata, x2i, i2x, NLFilename )

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2i, i2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    character(*), parameter :: subName = "(ice_init_mct) "
!-------------------------------------------------------------------------------


    if (present(NLFilename)) then
       call dice_comp_init(EClock, cdata, x2i, i2x, NLFilename)
    else
       call dice_comp_init(EClock, cdata, x2i, i2x)
    endif

end subroutine ice_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_run_mct
!
! !DESCRIPTION:
!     run method for dead ice model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ice_run_mct( EClock, cdata,  x2i, i2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2i        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: i2x        ! dead   -> driver

!EOP

   character(*), parameter :: subName = "(ice_run_mct) "
!-------------------------------------------------------------------------------

   call dice_comp_run(EClock, cdata, x2i, i2x)

end subroutine ice_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_final_mct
!
! !DESCRIPTION:
!     finalize method for dead ice model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine ice_final_mct( EClock, cdata, x2d, d2x)

   implicit none

   !----- arguments -----
   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d
   type(mct_aVect)             ,intent(inout) :: d2x

!EOP

   !--- formats ---
   character(*), parameter :: subName = "(ice_final_mct) "
!-------------------------------------------------------------------------------

   call dice_comp_final()

end subroutine ice_final_mct
!===============================================================================
!===============================================================================


end module ice_comp_mct
