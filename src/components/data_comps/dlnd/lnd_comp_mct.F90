module lnd_comp_mct

! !USES:

  use seq_cdata_mod
  use esmf
  use mct_mod

  use dlnd_comp_mod

! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_init_mct
!
! !DESCRIPTION:
!     initialize data lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_init_mct( EClock, cdata_l, x2l, l2x, NLFilename )

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_l
    type(mct_aVect)             , intent(inout) :: x2l, l2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    character(*), parameter :: subName = "(lnd_init_mct) "
!-------------------------------------------------------------------------------


    if (present(NLFilename)) then
       call dlnd_comp_init( EClock, cdata_l, x2l, l2x, NLFilename )
    else
       call dlnd_comp_init( EClock, cdata_l, x2l, l2x )
    endif

end subroutine lnd_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_run_mct
!
! !DESCRIPTION:
!     run method for data lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_run_mct( EClock, cdata_l,  x2l, l2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata_l
   type(mct_aVect)             ,intent(inout) :: x2l, l2x

!EOP

   character(*), parameter :: subName = "(lnd_run_mct) "
!-------------------------------------------------------------------------------

   call dlnd_comp_run( EClock, cdata_l, x2l, l2x)

end subroutine lnd_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_final_mct
!
! !DESCRIPTION:
!     finalize method for data lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine lnd_final_mct( EClock, cdata, x2l, l2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2l, l2x

!EOP

   !--- formats ---
   character(*), parameter :: subName = "(lnd_final_mct) "
!-------------------------------------------------------------------------------

   call dlnd_comp_final()

end subroutine lnd_final_mct
!===============================================================================
!===============================================================================


end module lnd_comp_mct
