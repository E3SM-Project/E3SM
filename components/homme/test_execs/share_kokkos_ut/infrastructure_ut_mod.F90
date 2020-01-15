module infrastructure_ut_mod

  implicit none

  public :: run_bfb_pow_f90

contains

  subroutine run_bfb_pow_f90 (x,e,y,n) bind(c)
    use iso_c_binding, only: c_double, c_int
    use bfb_mod, only: bfb_pow

    real (kind=c_double), intent(in)  :: x(n)
    real (kind=c_double), intent(in)  :: e
    real (kind=c_double), intent(out) :: y(n)
    integer (kind=c_int), intent(in)  :: n
    
    y = bfb_pow (x,e)
  end subroutine run_bfb_pow_f90
end module infrastructure_ut_mod
