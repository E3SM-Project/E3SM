module gw_iso_c
  use iso_c_binding
  implicit none

#include "eamxx_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to gw fortran.
!

contains

  subroutine gw_init_c(pver_in, pgwv_in, dc_in, cref_in, do_molec_diff_in, tau_0_ubc_in, nbot_molec_in, ktop_in, kbotbg_in, fcrit2_in, kwv_in, gravit_in, rair_in, alpha_in) bind(C)
    use gw_common, only : gw_common_init

    integer(kind=c_int) , value, intent(in) :: pver_in, pgwv_in, nbot_molec_in, ktop_in, kbotbg_in
    real(kind=c_real) , value, intent(in) :: dc_in, fcrit2_in, kwv_in, gravit_in, rair_in
    real(kind=c_real) , intent(in), dimension(-pgwv_in:) :: cref_in
    logical(kind=c_bool) , value, intent(in) :: do_molec_diff_in, tau_0_ubc_in
    real(kind=c_real) , intent(in), dimension(0:) :: alpha_in

    character(len=128) :: errstring

    call gw_common_init(pver_in, pgwv_in, dc_in, cref_in, do_molec_diff_in, tau_0_ubc_in, nbot_molec_in, ktop_in, kbotbg_in, fcrit2_in, kwv_in, gravit_in, rair_in, alpha_in, errstring)
  end subroutine gw_init_c

  subroutine gwd_compute_tendencies_from_stress_divergence_c(ncol, ngwv, do_taper, dt, effgw, tend_level, lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw) bind(C)
    use gw_common, only : gwd_compute_tendencies_from_stress_divergence, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: ncol, ngwv
    logical(kind=c_bool) , value, intent(in) :: do_taper
    real(kind=c_real) , value, intent(in) :: dt, effgw
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol) :: lat, xv, yv
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: dpm, rdpm, ubm, t, nm
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(inout), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver, -ngwv:ngwv) :: gwut
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: utgw, vtgw

    call gwd_compute_tendencies_from_stress_divergence(ncol, ngwv, do_taper, dt, effgw, tend_level, lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw)
  end subroutine gwd_compute_tendencies_from_stress_divergence_c
end module gw_iso_c
