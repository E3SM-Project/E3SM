!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_InitMod

!BOP
! !MODULE: POP_InitMod
! !DESCRIPTION:
!  This module contains the POP initialization method and initializes
!  everything needed by a POP simulation.  Primarily it is a driver
!  that calls individual initialization routines for each POP module.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id: POP_InitMod.F90 808 2006-04-28 17:06:38Z njn01 $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use initial

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_Initialize

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_Initialize
! !INTERFACE:

 subroutine POP_Initialize(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      errorCode              ! Returns an error code if any init fails

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  initialize return flag
!
!-----------------------------------------------------------------------

   ErrorCode = POP_Success

!-----------------------------------------------------------------------
!
!  call pop initialization routine
!
!-----------------------------------------------------------------------

   call initialize_POP

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize

!***********************************************************************

 end module POP_InitMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
