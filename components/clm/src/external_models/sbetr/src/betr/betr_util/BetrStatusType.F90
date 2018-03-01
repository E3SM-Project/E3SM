module BetrStatusType

use betr_constants, only :  betr_errmsg_len
implicit none

private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_status_type
    integer :: error
    character(len=betr_errmsg_len) :: msg
  contains
    procedure, public :: reset
    procedure, public :: set_msg
    procedure, public :: check_status
    procedure, public :: print_err
    procedure, public :: print_msg
  end type betr_status_type

  public :: create_betr_status_type
contains
!-------------------------------------------------------------------------------

  function create_betr_status_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_status_type), pointer :: create_betr_status_type
    class(betr_status_type), pointer :: betr_status

    allocate(betr_status)
    create_betr_status_type => betr_status

  end function create_betr_status_type
!-------------------------------------------------------------------------------
  subroutine reset(this)
  implicit none
  class(betr_status_type), intent(inout)  :: this

  this%msg = ''
  this%error = 0
  end subroutine reset
!-------------------------------------------------------------------------------
  subroutine set_msg(this, msg, err)
  implicit none
  class(betr_status_type), intent(inout)  :: this
  character(len=*), intent(in) :: msg
  integer, intent(in) :: err
  this%msg = trim(msg)
  this%error = err
  end subroutine set_msg
!-------------------------------------------------------------------------------
  function check_status(this)result(ans)
  implicit none
  class(betr_status_type), intent(in)  :: this
  logical :: ans
  if(this%error < 0)then
    ans = .true.
  else
    ans = .false.
  endif
  end function check_status
!-------------------------------------------------------------------------------
  function print_msg(this)
  implicit none
  class(betr_status_type), intent(in)  :: this
  character(len=betr_errmsg_len) :: print_msg

  print_msg = this%msg
  end function print_msg

!-------------------------------------------------------------------------------
  function print_err(this)
  implicit none
  class(betr_status_type), intent(in)  :: this
  integer :: print_err

  print_err = this%error
  end function print_err
end module BetrStatusType
