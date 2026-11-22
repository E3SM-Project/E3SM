module api_functions_mod
   use iso_c_binding, only: c_ptr, c_char
   implicit none
   public


contains
   subroutine to_c_string(str, cstr)
      character(len=*), intent(in) :: str
      character(kind=c_char), allocatable, target, intent(out) :: cstr(:)

      integer :: i, n

      n = len_trim(str)
      allocate (cstr(n + 1))

      do i = 1, n
         cstr(i) = str(i:i)
      end do
      cstr(n + 1) = c_null_char
   end subroutine to_c_string

   interface
      function c_harness_create() result(handle) bind(C, name="harness_create")
         import :: c_ptr
         type(c_ptr) :: handle
      end function

      subroutine c_harness_destroy(handle) bind(C, name="harness_destroy")
         import :: c_ptr
         type(c_ptr), value :: handle
      end subroutine

      subroutine c_harness_register_field(handle, name, extent, dtype, ptr) &
         bind(C, name="harness_register_field")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: handle, ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: extent, dtype
      end subroutine
   end interface

end module api_functions_mod
