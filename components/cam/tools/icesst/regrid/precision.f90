module precision

! Define 4-byte, 8-byte, and 16-byte sizes

  integer, parameter :: r4 = selected_real_kind(5)
  integer, parameter :: r8 = selected_real_kind(12)
  integer, parameter :: r16 = selected_real_kind(20)
  integer, parameter :: i8 = selected_int_kind(13)
end module precision
