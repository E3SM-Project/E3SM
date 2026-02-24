module emulator_handle_mod
  use, intrinsic :: iso_c_binding
  use emulator_f2c_api
  implicit none

  type :: emulator_handle
     type(c_ptr) :: h = c_null_ptr
  contains
     final :: destroy_emulator_handle
  end type emulator_handle

contains

  subroutine destroy_emulator_handle(self)
    type(emulator_handle), intent(inout) :: self
    if (c_associated(self%h)) then
       call emulator_finalize(self%h)
       self%h = c_null_ptr
    end if
  end subroutine destroy_emulator_handle

end module emulator_handle_mod
