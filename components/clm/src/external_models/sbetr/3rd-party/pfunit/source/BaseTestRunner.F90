!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: BaseTestRunner
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
module BaseTestRunner_mod
   use TestListener_mod
   implicit none
   private

   public :: BaseTestRunner

   type, abstract, extends(TestListener) :: BaseTestRunner
      private

   contains
      procedure(run2), deferred :: run
   end type BaseTestRunner

   abstract interface

      ! TODO - report bug to NAG.  If this is named "run" then
      ! RubustRunner fails to compile with message about conflicting types

      function run2(this, aTest, context) result(result)
         use Test_mod
         use ParallelContext_mod
         use TestResult_mod
         import BaseTestRunner

         type (TestResult) :: result
         class (BaseTestRunner), intent(inout) :: this
         class (Test), intent(inout) :: aTest
         class (ParallelContext), intent(in) :: context
      end function run2

   end interface

end module BaseTestRunner_mod
