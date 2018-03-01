module prec
!
! Define 8-byte size
! Use 6 for real*4 arithmetic.  Reason is to get bfb agreement with
! Karl Taylor's original code.
!
!  integer, parameter :: r8 = selected_real_kind(6)
  integer, parameter :: r8 = selected_real_kind(12)
end module prec
