module ekat_string_utils

  implicit none

contains

  subroutine string_c2f (c_str_ptr, f_str)
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR, c_ptr, c_f_pointer
    !
    ! Input(s)
    !
    type(c_ptr),        intent(in)  :: c_str_ptr
    character(len=256), intent(out) :: f_str
    !
    ! Local(s)
    !
    character(kind=C_CHAR,len=256), pointer :: c_str 
    integer :: str_len

    call c_f_pointer(c_str_ptr,c_str)
    str_len = index(c_str, C_NULL_CHAR) - 1

    f_str = trim(c_str(1:str_len))
  end subroutine string_c2f

  subroutine string_f2c (f_str, c_str)
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR
    !
    ! Input(s)
    !
    character(len=*),               intent(in)  :: f_str
    character(kind=C_CHAR,len=256), intent(out) :: c_str

    c_str = TRIM(f_str)//C_NULL_CHAR
  end subroutine string_f2c

end module ekat_string_utils
