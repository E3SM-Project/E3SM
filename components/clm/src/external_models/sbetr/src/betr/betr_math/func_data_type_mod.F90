module func_data_type_mod
  ! NOTE(bja, 201604) this type needs to be defined in a separate
  ! module to avoid a dependency problems.

  use bshr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! Kitchen sink container that can be used to pass data into call
  ! back functions. To keep the size down:
  !
  ! * all large arrays should be dynamically allocated.
  !
  ! * variables in the struct should be kept generic, then use
  ! assaciate blocks in the calling routine and call back functions to
  ! give them meaningful names.
  !
  type, public:: func_data_type
     real(r8), pointer :: aj(:) => null()
     real(r8) :: iJ
     integer :: nJ
  end type func_data_type

end module func_data_type_mod
