!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: BrokenSetUpCase_mod
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
module BrokenSetUpCase_mod
   use TestCase_mod, only: TestCase
   use Exception_mod, only: throw
   implicit none
   private
   
   public :: BrokenSetUpCase
   public :: newBrokenSetUpCase
   
   type, extends(TestCase) :: BrokenSetUpCase
      private
      character(len=40), public :: runLog
   contains
      procedure :: setUp
      procedure :: runMethod
      procedure :: tearDown
   end type BrokenSetUpCase
   
contains

   function newBrokenSetUpCase() result(this)
      type (BrokenSetUpCase), pointer :: this
      allocate(this)
      call this%setName('BrokenSetUpCase')
   end function newBrokenSetUpCase

   subroutine setUp(this)
      class(BrokenSetUpCase), intent(inOut) :: this

      this%runLog = 'broken setUp'
      call throw('This setUp() is intentionally broken.')

   end subroutine setUp

   subroutine tearDown(this)
      class(BrokenSetUpCase), intent(inOut) :: this

      this%runLog = trim(this%runLog)//' tearDown'

   end subroutine tearDown

   subroutine runMethod(this)
      class(BrokenSetUpCase), intent(inOut) :: this

      this%runLog = trim(this%runLog)//' run'

   end subroutine runMethod

end module BrokenSetUpCase_mod
