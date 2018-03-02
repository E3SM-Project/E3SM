!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: FixtureTestCase_mod
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
module FixtureTestCase_mod
   use TestCase_mod, only: TestCase
   implicit none
   private

   public :: FixtureTestCase
   public :: newFixtureTestCase
   public :: delete
   public :: methodA
   public :: methodB
   public :: simpleTestMethod

   type, extends(TestCase) :: FixtureTestCase
      private
      character(len=30), public :: runLog
   contains
      procedure :: setUp
      procedure :: runMethod
      procedure :: tearDown
   end type FixtureTestCase

   interface delete
      module procedure delete_
   end interface
   
contains

   function newFixtureTestCase() result(this)
      type(FixtureTestCase) :: this

      call this%setName('FixtureTestCase')
      this%runLog = ' '

   end function newFixtureTestCase

   subroutine setUp(this)
      class(FixtureTestCase), intent(inOut) :: this
      this%runLog = trim(this%runLog) // 'setUp '
   end subroutine setUp

   subroutine tearDown(this)
      class(FixtureTestCase), intent(inOut) :: this
      this%runLog = trim(this%runLog) // ' tearDown'
   end subroutine tearDown

   subroutine runMethod(this)
      class(FixtureTestCase), intent(inOut) :: this

      this%runLog = trim(this%runLog) // ' run'

   end subroutine runMethod

   subroutine simpleTestMethod(this)
      class (FixtureTestCase), intent(inOut) :: this
      this%runLog = trim(this%runLog) // ' run'
   end subroutine simpleTestMethod

   subroutine methodA(this)
      class (FixtureTestCase), intent(inOut) :: this
      this%runLog = trim(this%runLog) // ' methodA'
   end subroutine methodA

   subroutine methodB(this)
      class (FixtureTestCase), intent(inOut) :: this
      this%runLog = trim(this%runLog) // ' methodB'
   end subroutine methodB

   subroutine delete_(this)
      type (FixtureTestCase), intent(inOut) :: this
   end subroutine delete_

end module FixtureTestCase_mod
