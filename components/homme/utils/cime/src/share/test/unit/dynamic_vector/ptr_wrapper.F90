module ptr_wrapper

! This module simply defines a wrapper for integer pointers in order to
! test the dynamic_vector type on a derived type.

implicit none
private
save

public :: int_ptr

type int_ptr
   integer, pointer :: p => null()
 contains
   procedure, pass(first) :: cmp => int_ptr_cmp
   generic :: operator(==) => cmp
end type int_ptr

contains

elemental function int_ptr_cmp(first, second) result(is_same)
  class(int_ptr), intent(in) :: first
  class(int_ptr), intent(in) :: second
  logical :: is_same

  is_same = associated(first%p, second%p) .or. &
       (.not. associated(first%p) .and. .not. associated(second%p))

end function int_ptr_cmp

end module ptr_wrapper
