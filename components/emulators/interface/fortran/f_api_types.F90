module f_api_types_mod
  implicit none
  private

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


end module f_api_types_mod
