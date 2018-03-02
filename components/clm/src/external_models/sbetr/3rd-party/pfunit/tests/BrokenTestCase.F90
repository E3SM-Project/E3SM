!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: BrokenTestCase_mod
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune,  NASA/GSFC
!!
!! @date
!! 20 Mar 2015
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 20 Mar 2015 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------
module BrokenTestCase_mod
   use TestCase_mod, only: TestCase
   implicit none
   private
   
   public :: BrokenTestCase
   
   type, extends(TestCase) :: BrokenTestCase
      private
      character(len=40), public :: runLog
   contains
      procedure :: setUp
      procedure :: tearDown
      procedure :: runMethod
   end type BrokenTestCase
   
contains

   subroutine setUp(this)
      class(BrokenTestCase), intent(inOut) :: this

      this%runLog = 'setUp'

   end subroutine setUp

   subroutine runMethod(this)
      use Exception_mod, only: throw
      class(BrokenTestCase), intent(inOut) :: this

      this%runLog = trim(this%runLog) // ' broken run'
      call throw('This test is intentionally broken.')

   end subroutine runMethod

   subroutine tearDown(this)
      class(BrokenTestCase), intent(inOut) :: this

      this%runLog = trim(this%runLog) // ' tearDown'

   end subroutine tearDown

end module BrokenTestCase_mod
