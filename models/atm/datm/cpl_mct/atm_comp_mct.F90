module atm_comp_mct

! !USES:

  use seq_cdata_mod
  use esmf
  use mct_mod

  use datm_comp_mod

! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_init_mct
!
! !DESCRIPTION:
!     initialize data atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_init_mct( EClock, cdata, x2a, a2x, NLFilename )

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2a, a2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    character(*), parameter :: subName = "(atm_init_mct) "
!-------------------------------------------------------------------------------


    if (present(NLFilename)) then
       call datm_comp_init(EClock, cdata, x2a, a2x, NLFilename)
    else
       call datm_comp_init(EClock, cdata, x2a, a2x)
    endif

end subroutine atm_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_run_mct
!
! !DESCRIPTION:
!     run method for dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_run_mct( EClock, cdata,  x2a, a2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2a        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: a2x        ! dead   -> driver

!EOP

   character(*), parameter :: subName = "(atm_run_mct) "
!-------------------------------------------------------------------------------

   call datm_comp_run(EClock, cdata, x2a, a2x)

end subroutine atm_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_final_mct
!
! !DESCRIPTION:
!     finalize method for dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine atm_final_mct(EClock, cdata, x2d, d2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

!EOP

   !--- formats ---
   character(*), parameter :: subName = "(atm_final_mct) "
!-------------------------------------------------------------------------------

   call datm_comp_final()

end subroutine atm_final_mct
!===============================================================================
!===============================================================================


end module atm_comp_mct
