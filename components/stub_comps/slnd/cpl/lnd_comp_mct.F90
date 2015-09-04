module lnd_comp_mct

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

  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_init_mct
!
! !DESCRIPTION:
!     stub lnd model init
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine lnd_init_mct( EClock, cdata, x2d, d2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename 

!EOP
!-------------------------------------------------------------------------------

   call seq_infodata_PutData(cdata%infodata, &
        lnd_present=.false., lnd_prognostic=.false.)

end subroutine lnd_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_run_mct
!
! !DESCRIPTION:
!     stub lnd model run
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine lnd_run_mct( EClock, cdata, x2d, d2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        
   type(mct_aVect)             ,intent(inout) :: d2x        

!EOP
!-------------------------------------------------------------------------------

end subroutine lnd_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: lnd_final_mct
!
! !DESCRIPTION:
!     stub lnd model finalize 
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine lnd_final_mct( EClock, cdata, x2d, d2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        
    type(mct_aVect)             ,intent(inout) :: d2x        

!EOP
!-------------------------------------------------------------------------------

 end subroutine lnd_final_mct

!===============================================================================

end module lnd_comp_mct
