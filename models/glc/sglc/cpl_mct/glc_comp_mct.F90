module glc_comp_mct

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

  public :: glc_init_mct
  public :: glc_run_mct
  public :: glc_final_mct

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_init_mct
!
! !DESCRIPTION:
!     stub glc model init
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine glc_init_mct( EClock, cdata, x2d, d2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename

!EOP
!-------------------------------------------------------------------------------

   call seq_infodata_PutData(cdata%infodata, glc_present=.false., &
        glc_prognostic=.false.)

end subroutine glc_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_run_mct
!
! !DESCRIPTION:
!     stub glc model run
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine glc_run_mct( EClock, cdata, x2d, d2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        
   type(mct_aVect)             ,intent(inout) :: d2x        

!EOP
!-------------------------------------------------------------------------------

end subroutine glc_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_final_mct
!
! !DESCRIPTION:
!     stub glc model finalize 
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine glc_final_mct( EClock, cdata, x2d, d2x)

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        
   type(mct_aVect)             ,intent(inout) :: d2x        

!EOP
!-------------------------------------------------------------------------------

 end subroutine glc_final_mct

!===============================================================================

end module glc_comp_mct
