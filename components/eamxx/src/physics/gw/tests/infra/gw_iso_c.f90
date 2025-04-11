module gw_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to gw fortran.
!

contains

  subroutine gwd_compute_tendencies_from_stress_divergence_c(ncol, pver, pgwv, ngwv, do_taper, dt, effgw, tend_level, lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw) bind(C)
    use gw, only : gwd_compute_tendencies_from_stress_divergence

    integer(kind=c_int) , value, intent(in) :: ncol, pver, pgwv, ngwv
    logical(kind=c_bool) , value, intent(in) :: do_taper
    real(kind=c_real) , value, intent(in) :: dt, effgw
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol) :: lat, xv, yv
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: dpm, rdpm, ubm, t, nm
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(inout), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver, -ngwv:ngwv) :: gwut
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: utgw, vtgw

    call gwd_compute_tendencies_from_stress_divergence(ncol, pver, pgwv, ngwv, do_taper, dt, effgw, tend_level, lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw)
  end subroutine gwd_compute_tendencies_from_stress_divergence_c
end module gw_iso_c
