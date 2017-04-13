module rof_comp_mct

! !USES:

  use seq_cdata_mod
  use esmf
  use mct_mod

  use drof_comp_mod

! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: rof_init_mct
  public :: rof_run_mct
  public :: rof_final_mct

!--------------------------------------------------------------------------
! Private data
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
!     initialize data rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine rof_init_mct( EClock, cdata, x2r, r2x, NLFilename )

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2r, r2x
    character(len=*), optional  , intent(in)    :: NLFilename

!EOP

    character(*), parameter :: subName = "(rof_init_mct) "
!-------------------------------------------------------------------------------


    if (present(NLFilename)) then
       call drof_comp_init( EClock, cdata, x2r, r2x, NLFilename )
    else
       call drof_comp_init( EClock, cdata, x2r, r2x)
    endif

end subroutine rof_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_run_mct
!
! !DESCRIPTION:
!     run method for data rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine rof_run_mct( EClock, cdata, x2r, r2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2r
   type(mct_aVect)             ,intent(inout) :: r2x

!EOP

   character(*), parameter :: subName = "(rof_run_mct) "
!-------------------------------------------------------------------------------

   call drof_comp_run( EClock, cdata, x2r, r2x)

end subroutine rof_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_final_mct
!
! !DESCRIPTION:
!     finalize method for data rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine rof_final_mct( EClock, cdata, x2r, r2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2r
    type(mct_aVect)             ,intent(inout) :: r2x

!EOP

   !--- formats ---
   character(*), parameter :: subName = "(rof_final_mct) "
!-------------------------------------------------------------------------------

   call drof_comp_final()

end subroutine rof_final_mct
!===============================================================================
!===============================================================================


end module rof_comp_mct
