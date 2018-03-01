module BetrStatusSimType
use BeTRStatusType, only : betr_status_type
implicit none

private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

type, public, extends(betr_status_type) :: betr_status_sim_type

  integer :: cindex
contains
  procedure, public :: setcol

end type betr_status_sim_type
 public :: create_betr_status_sim_type
contains

!-------------------------------------------------------------------------------

  function create_betr_status_sim_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_status_sim_type), pointer :: create_betr_status_sim_type
    class(betr_status_sim_type), pointer :: betr_status

    allocate(betr_status)
    create_betr_status_sim_type => betr_status

  end function create_betr_status_sim_type

!-------------------------------------------------------------------------------

  subroutine setcol(this, cc)
  implicit none
  class(betr_status_sim_type) :: this
  integer :: cc
  this%cindex = cc

  end subroutine setcol

!-------------------------------------------------------------------------------

end module BetrStatusSimType
