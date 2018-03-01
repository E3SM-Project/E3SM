!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Assert
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune, NASA/GSFC 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Nov 2013 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------
module Assert_mod
   use AssertBasic_mod
#include "AssertArrays.fh"
   implicit none
   private

   public :: assertFail
   public :: assertTrue
   public :: assertFalse
   public :: assertEqual
   public :: assertExceptionRaised
   public :: assertSameShape

   public :: assertAny
   public :: assertAll
   public :: assertNone
   public :: assertNotAll

   public :: assertNotEqual
   public :: assertLessThan, assertLessThanOrEqual
   public :: assertGreaterThan, assertGreaterThanOrEqual
   public :: assertRelativelyEqual

   public :: assertIsNan, assertIsFinite

   ! Optional arguments for assertEqual.
   public :: WhitespaceOptions
   public :: IGNORE_ALL, TRIM_ALL, KEEP_ALL, IGNORE_DIFFERENCES

contains

end module Assert_mod
