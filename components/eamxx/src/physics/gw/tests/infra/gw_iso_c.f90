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

  subroutine gw_init_c(pver_in, pgwv_in, dc_in, cref_in, orographic_only, do_molec_diff_in, tau_0_ubc_in, nbot_molec_in, ktop_in, kbotbg_in, fcrit2_in, kwv_in, gravit_in, rair_in, alpha_in) bind(C)
    use gw_common, only : gw_common_init

    integer(kind=c_int) , value, intent(in) :: pver_in, pgwv_in, nbot_molec_in, ktop_in, kbotbg_in
    real(kind=c_real) , value, intent(in) :: dc_in, fcrit2_in, kwv_in, gravit_in, rair_in
    real(kind=c_real) , intent(in), dimension(-pgwv_in:pgwv_in) :: cref_in
    logical(kind=c_bool) , value, intent(in) :: orographic_only, do_molec_diff_in, tau_0_ubc_in
    real(kind=c_real) , intent(in), dimension(0:pver_in) :: alpha_in

    character(len=128) :: errstring

    call gw_common_init(pver_in, pgwv_in, dc_in, cref_in, orographic_only, do_molec_diff_in, tau_0_ubc_in, nbot_molec_in, ktop_in, kbotbg_in, fcrit2_in, kwv_in, gravit_in, rair_in, alpha_in, errstring)
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

  subroutine gw_prof_c(ncol, cpair, t, pmid, pint, rhoi, ti, nm, ni) bind(C)
    use gw_common, only : gw_prof, pver

    integer(kind=c_int) , value, intent(in) :: ncol
    real(kind=c_real) , value, intent(in) :: cpair
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: t, pmid
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: pint
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: rhoi, ti, ni
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: nm

    call gw_prof(ncol, cpair, t, pmid, pint, rhoi, ti, nm, ni)
  end subroutine gw_prof_c

  subroutine momentum_energy_conservation_c(ncol, tend_level, dt, taucd, pint, pdel, u, v, dudt, dvdt, dsdt, utgw, vtgw, ttgw) bind(C)
    use gw_common, only : momentum_energy_conservation, pver

    integer(kind=c_int) , value, intent(in) :: ncol
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , value, intent(in) :: dt
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver, 4) :: taucd
    real(kind=c_real) , intent(in), dimension(ncol, pver+1) :: pint
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: pdel, u, v
    real(kind=c_real) , intent(inout), dimension(ncol, pver) :: dudt, dvdt, dsdt
    real(kind=c_real) , intent(inout), dimension(ncol, pver) :: utgw, vtgw, ttgw

    call momentum_energy_conservation(ncol, tend_level, dt, taucd, pint, pdel, u, v, dudt, dvdt, dsdt, utgw, vtgw, ttgw)
  end subroutine momentum_energy_conservation_c

  subroutine gwd_compute_stress_profiles_and_diffusivities_c(ncol, ngwv, src_level, ubi, c, rhoi, ni, kvtt, t, ti, piln, tau) bind(C)
    use gw_common, only : gwd_compute_stress_profiles_and_diffusivities, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: ncol, ngwv
    integer(kind=c_int) , intent(in), dimension(ncol) :: src_level
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: ubi, rhoi, ni, kvtt, ti, piln
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: t
    real(kind=c_real) , intent(inout), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau

    call gwd_compute_stress_profiles_and_diffusivities(ncol, ngwv, src_level, ubi, c, rhoi, ni, kvtt, t, ti, piln, tau)
  end subroutine gwd_compute_stress_profiles_and_diffusivities_c

  subroutine gwd_project_tau_c(ncol, ngwv, tend_level, tau, ubi, c, xv, yv, taucd) bind(C)
    use gw_common, only : gwd_project_tau, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: ncol, ngwv
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: ubi
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(in), dimension(ncol) :: xv, yv
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver, 4) :: taucd

    call gwd_project_tau(ncol, ngwv, tend_level, tau, ubi, c, xv, yv, taucd)
  end subroutine gwd_project_tau_c
end module gw_iso_c
