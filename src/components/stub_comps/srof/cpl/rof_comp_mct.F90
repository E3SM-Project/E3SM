module rof_comp_mct

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

  public :: rof_init_mct
  public :: rof_run_mct
  public :: rof_final_mct
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_init_mct
!
! !DESCRIPTION:
!     stub rof model init
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine rof_init_mct( EClock, cdata, x2r, r2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2r, r2x
    character(len=*), optional  , intent(in)    :: NLFilename

!EOP
!-------------------------------------------------------------------------------

   call seq_infodata_PutData(cdata%infodata, rof_present=.false., &
      rofice_present=.false., rof_prognostic=.false.)
   call seq_infodata_PutData(cdata%infodata, flood_present=.false.)

end subroutine rof_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_run_mct
!
! !DESCRIPTION:
!     stub rof model run
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine rof_run_mct( EClock, cdata, x2r, r2x )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2r, r2x

!EOP
!-------------------------------------------------------------------------------

end subroutine rof_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: rof_final_mct
!
! !DESCRIPTION:
!     stub rof model finalize
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine rof_final_mct( EClock, cdata, x2r, r2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2r, r2x

!EOP
!-------------------------------------------------------------------------------

 end subroutine rof_final_mct

!===============================================================================

end module rof_comp_mct
