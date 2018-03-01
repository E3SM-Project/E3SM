!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: MpiTestMethod
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
module MpiTestMethod_mod
   use Test_mod
   use TestCase_mod
   use MpiTestCase_mod
   use MpiTestParameter_mod
   implicit none
   private

   public :: MpiTestMethod
   public :: newMpiTestMethod

   interface newMpiTestMethod
      module procedure newMpiTest_basic
      module procedure newMpiTest_setUpTearDown
   end interface newMpiTestMethod

!!$   interface MpiTestMethod
!!$      module procedure newMpiTest_basic
!!$      module procedure newMpiTest_setUpTearDown
!!$   end interface MpiTestMethod

   type, extends(MpiTestCase) :: MpiTestMethod
      procedure(mpiMethod), pointer :: userMethod => null()
      procedure(mpiMethod), nopass, pointer :: userSetUp => null()
      procedure(mpiMethod), nopass, pointer :: userTearDown => null()
   contains
      procedure :: runMethod
      procedure :: setUp
      procedure :: tearDown
   end type MpiTestMethod

   abstract interface
      subroutine mpiMethod(this)
         import MpiTestMethod
         class (MpiTestMethod), intent(inout) :: this
      end subroutine mpiMethod
   end interface

contains

   function newMpiTest_basic(name, userMethod, numProcesses) result(mpiTest)
      character(len=*), intent(in) :: name
      procedure (runMethod) :: userMethod
      integer, intent(in) :: numProcesses
      type (MpiTestMethod) :: mpiTest

      call mpiTest%setName(name)
      mpiTest%userMethod => userMethod
      call mpiTest%setTestParameter(MpiTestParameter(numProcesses))

   end function newMpiTest_basic

   function newMpiTest_setUpTearDown(name, userMethod, numProcesses, setUp, tearDown) result(mpiTest)
      character(len=*), intent(in) :: name
      procedure (runMethod) :: userMethod
      integer, intent(in) :: numProcesses
      type (MpiTestMethod) :: mpiTest
      procedure (runMethod) :: setUp
      procedure (runMethod) :: tearDown

      call mpiTest%setName(name)
      mpiTest%userMethod => userMethod
      call mpiTest%setTestParameter(MpiTestParameter(numProcesses))

      mpiTest%userSetUp => setUp
      mpiTest%userTearDown => tearDown

   end function newMpiTest_setUpTearDown

   subroutine runMethod(this)
      class (MpiTestMethod), intent(inout) :: this
      call this%userMethod()
   end subroutine runMethod

   subroutine setUp(this)
      class (MpiTestMethod), intent(inout) :: this
      if (associated(this%userSetUp)) call this%userSetUp(this)
   end subroutine setUp

   subroutine tearDown(this)
      class (MpiTestMethod), intent(inout) :: this
      if (associated(this%userTearDown)) call this%userTearDown(this)
   end subroutine tearDown

end module MpiTestMethod_mod
