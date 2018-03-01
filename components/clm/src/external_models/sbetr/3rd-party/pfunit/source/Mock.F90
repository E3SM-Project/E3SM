
!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Mock
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune, NASA/GSFC
!!
!! @date
!! 12 May 2014
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 12 May 2014 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------

module Mock_mod

  use MockRepository_mod

  implicit none
  private

  public :: Mock

  type Mock
     character(len=MAX_LEN_METHOD_NAME) :: name = ''
     type(MockRepository) :: repository
!     procedure(subI), pointer, nopass :: sub
   contains
     procedure registerCall
  end type Mock


contains

!  function newMock(name,signature,mockRepo) result(mock)

  subroutine registerCall(this)
    class(Mock), intent(inout) :: this
    call this%repository%registerMockCallBy(this%name)
  end subroutine registerCall

end module Mock_mod


