!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: ParameterizedTestCase
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
! 01 Jan 2014 - Added "hidden" method toStringActual()
!
!-------------------------------------------------------------------------------
module ParameterizedTestCase_mod
   use TestCase_mod
   use AbstractTestParameter_mod
   implicit none
   private
   
   public :: ParameterizedTestCase
   public :: MAX_LEN_LABEL

   integer, parameter :: MAX_LEN_LABEL = 32
   type, abstract, extends(TestCase) :: ParameterizedTestCase
      class (AbstractTestParameter), allocatable :: testParameter
   contains
      procedure :: getName ! override from TestCase
      procedure :: setTestParameter
   end type ParameterizedTestCase

contains


   function getName(this) result(name)
      class (ParameterizedTestCase), intent(in) :: this
      character(:), allocatable :: name

      name = trim(this%baseName()) // '[' // trim(this%testParameter%toStringActual()) // ']'

   end function getName

   subroutine setTestParameter(this, testParameter)
      class (ParameterizedTestCase), intent(inout) :: this
      class (AbstractTestParameter), intent(in) :: testParameter
      allocate(this%testParameter, source=testParameter)
   end subroutine setTestParameter

end module ParameterizedTestCase_mod
