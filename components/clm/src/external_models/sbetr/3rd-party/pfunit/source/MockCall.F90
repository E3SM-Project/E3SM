!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: MockCall
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
module MockCall_mod
   use Exception_mod
   implicit none
   private

   public :: MockCall
   public :: newMockCall
   integer, parameter :: MAXLEN_METHOD_NAME = 32
   type MockCall
      character(len=MAXLEN_METHOD_NAME) :: methodName
      class(*), pointer :: argument
   contains
      procedure :: expect
      procedure :: getExpectedValue
   end type MockCall

contains

   function newMockCall(name) result(mCall)
      character(len=*), intent(in) :: name
      type (MockCall) :: mCall

      mCall%methodName = name
   end function NewMockCall

   subroutine expect(this, expectedArgument)
      class (MockCall), intent(inout) :: this
      class(*), target, intent(in) :: expectedArgument
      this%argument => expectedArgument
   end subroutine expect

   function getExpectedValue(this) result(p)
      class(MockCall), intent(in) :: this
      class(*), pointer :: p

      p => this%argument
   end function getExpectedValue

end module MockCall_mod
