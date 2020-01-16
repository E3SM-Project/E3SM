module infrastructure_ut_mod

  implicit none

  public :: power_of_two_f90, int_pow_f90, bfb_pow_f90

contains

  subroutine power_of_two_f90 (e,y) bind(c)
    use iso_c_binding, only: c_double
    use bfb_mod, only: power_of_two

    real (kind=c_double), intent(in)  :: e
    real (kind=c_double), intent(out) :: y
    
    y = power_of_two (e)
  end subroutine power_of_two_f90

  subroutine int_pow_f90 (x,e,y) bind(c)
    use iso_c_binding, only: c_double, c_int
    use bfb_mod, only: int_pow

    real (kind=c_double), intent(in)  :: x
    integer (kind=c_int), intent(in)  :: e
    real (kind=c_double), intent(out) :: y
    
    y = int_pow (x,e)
  end subroutine int_pow_f90

  subroutine bfb_pow_f90 (x,e,y,n) bind(c)
    use iso_c_binding, only: c_double, c_int
    use bfb_mod, only: bfb_pow

    real (kind=c_double), intent(in)  :: x(n)
    real (kind=c_double), intent(in)  :: e
    real (kind=c_double), intent(out) :: y(n)
    integer (kind=c_int), intent(in)  :: n
    
    y = bfb_pow (x,e)
  end subroutine bfb_pow_f90
end module infrastructure_ut_mod
