module tms_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from tms fortran to scream c++.
!

interface

  subroutine compute_tms_f(ncols, nlevs, u_wind, v_wind, t_mid, p_mid, &
                           exner, z_mid, sgh, landfrac, ksrf, taux, tauy) bind (C)
    use iso_c_binding

    integer(kind=c_int), intent(in), value :: ncols, nlevs
    real(kind=c_real), intent(in) :: u_wind(ncols,nlevs)
    real(kind=c_real), intent(in) :: v_wind(ncols,nlevs)
    real(kind=c_real), intent(in) :: t_mid(ncols,nlevs)
    real(kind=c_real), intent(in) :: p_mid(ncols,nlevs)
    real(kind=c_real), intent(in) :: exner(ncols,nlevs)
    real(kind=c_real), intent(in) :: z_mid(ncols,nlevs)

    real(kind=c_real), intent(in) :: sgh(ncols)
    real(kind=c_real), intent(in) :: landfrac(ncols)

    real(kind=c_real), intent(out) :: ksrf(ncols)
    real(kind=c_real), intent(out) :: taux(ncols)
    real(kind=c_real), intent(out) :: tauy(ncols)

  end subroutine compute_tms_f
end interface

end module tms_iso_f
