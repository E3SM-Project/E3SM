!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: SurrogateTestCase
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
!
! This module exists as part of the Surrogate design pattern which
! helps to circumvent circular dependenciens betwenn Fortran classes.
! In this case, the Test hierarchy depends upon TestResult.  In turn
! TestResult depends on TestCase, which is a subclass of Test.
!
! This is a modified variant of the Surrogate pattern due to the
! injection within an inheritance hierarchy (i.e. between Test and
! TestCase).  Since Fortran only supports single inheritance, the
! Multiple-Inheritance design pattern is also required.  That portion
! is implemented in the TestCase module.

module SurrogateTestCase_mod
   implicit none
   private

   public :: SurrogateTestCase

   type, abstract :: SurrogateTestCase
      private
   contains
      procedure(getName), deferred :: getName 
      procedure(setName), deferred :: setName
      procedure(runBare), deferred :: runBare
   end type SurrogateTestCase

   abstract interface

      ! Run the SUT and assert the results
      subroutine runBare(this)
         import SurrogateTestCase
         class (SurrogateTestCase), intent(inout) :: this
      end subroutine runBare
      
      ! Return the name for TestCase (may need to move to Test)
      function getName(this) result(name)
         import SurrogateTestCase
         class (SurrogateTestCase), intent(in) :: this
         character(:), allocatable :: name
      end function getName
      
      ! Set the test name for TestCase (may need to move to Test)
      subroutine setName(this, name)
         import SurrogateTestCase
         class (SurrogateTestCase), intent(inout) :: this
         character(len=*),intent(in) :: name
      end subroutine setName
      
   end interface

end module SurrogateTestCase_mod
