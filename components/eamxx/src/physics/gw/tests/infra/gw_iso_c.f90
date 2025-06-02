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

  subroutine gwd_precalc_rhoi_c(pcnst, ncol, ngwv, dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm, c, q, dse, egwdffi, qtgw, dttdf, dttke, ttgw) bind(C)
    use gw_common, only : gwd_precalc_rhoi, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: pcnst, ncol, ngwv
    real(kind=c_real) , value, intent(in) :: dt
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: pmid, t, ubm, nm, rdpm, dse
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: pint
    real(kind=c_real) , intent(in), dimension(ncol, pver, -ngwv:ngwv) :: gwut
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(in), dimension(ncol, pver, pcnst) :: q
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: egwdffi
    real(kind=c_real) , intent(out), dimension(ncol, pver, pcnst) :: qtgw
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: dttdf, dttke, ttgw

    call gwd_precalc_rhoi(ncol, ngwv, dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm, c, q, dse, egwdffi, qtgw, dttdf, dttke, ttgw)
  end subroutine gwd_precalc_rhoi_c

  subroutine gw_drag_prof_c(pcnst, ncol, ngwv, src_level, tend_level, do_taper, dt, lat, t, ti, pmid, pint, dpm, rdpm, piln, rhoi, nm, ni, ubm, ubi, xv, yv, effgw, c, kvtt, q, dse, tau, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, gwut, dttdf, dttke) bind(C)
    use gw_common, only : gw_drag_prof, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: pcnst, ncol, ngwv
    integer(kind=c_int) , intent(in), dimension(ncol) :: src_level, tend_level
    logical(kind=c_bool) , value, intent(in) :: do_taper
    real(kind=c_real) , value, intent(in) :: dt, effgw
    real(kind=c_real) , intent(in), dimension(ncol) :: lat, xv, yv
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: t, pmid, dpm, rdpm, nm, ubm, dse
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: ti, pint, piln, rhoi, ni, ubi, kvtt
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(in), dimension(ncol, pver, pcnst) :: q
    real(kind=c_real) , intent(inout), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: utgw, vtgw, ttgw, dttdf, dttke
    real(kind=c_real) , intent(out), dimension(ncol, pver, pcnst) :: qtgw
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver, 4) :: taucd
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: egwdffi
    real(kind=c_real) , intent(out), dimension(ncol, pver, -ngwv:ngwv) :: gwut

    call gw_drag_prof(ncol, ngwv, src_level, tend_level, do_taper, dt, lat, t, ti, pmid, pint, dpm, rdpm, piln, rhoi, nm, ni, ubm, ubi, xv, yv, effgw, c, kvtt, q, dse, tau, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, gwut, dttdf, dttke)
  end subroutine gw_drag_prof_c

  subroutine gw_front_init_c(taubgnd, frontgfc_in, kfront_in) bind(C)
    use gw_front, only : gw_front_init

    real(kind=c_real) , value, intent(in) :: taubgnd, frontgfc_in
    integer(kind=c_int) , value, intent(in) :: kfront_in

    character(len=128) :: errstring

    call gw_front_init(taubgnd, frontgfc_in, kfront_in, errstring)
  end subroutine gw_front_init_c

  subroutine gw_front_project_winds_c(ncol, kbot, u, v, xv, yv, ubm, ubi) bind(C)
    use gw_common, only : pver
    use gw_front, only : gw_front_project_winds

    integer(kind=c_int) , value, intent(in) :: ncol, kbot
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: u, v
    real(kind=c_real) , intent(out), dimension(ncol) :: xv, yv
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: ubm
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: ubi

    call gw_front_project_winds(ncol, kbot, u, v, xv, yv, ubm, ubi)
  end subroutine gw_front_project_winds_c

  subroutine gw_front_gw_sources_c(ncol, ngwv, kbot, frontgf, tau) bind(C)
    use gw_common, only : pver, pgwv
    use gw_front, only : gw_front_gw_sources

    integer(kind=c_int) , value, intent(in) :: ncol, ngwv, kbot
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: frontgf
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau

    call gw_front_gw_sources(ncol, ngwv, kbot, frontgf, tau)
  end subroutine gw_front_gw_sources_c

  subroutine gw_cm_src_c(ncol, ngwv, kbot, u, v, frontgf, src_level, tend_level, tau, ubm, ubi, xv, yv, c) bind(C)
    use gw_common, only : pver, pgwv
    use gw_front, only : gw_cm_src

    integer(kind=c_int) , value, intent(in) :: ncol, ngwv, kbot
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: u, v
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: frontgf
    integer(kind=c_int) , intent(out), dimension(ncol) :: src_level, tend_level
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: ubm
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: ubi
    real(kind=c_real) , intent(out), dimension(ncol) :: xv, yv
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv) :: c

    call gw_cm_src(ncol, ngwv, kbot, u, v, frontgf, src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  end subroutine gw_cm_src_c
end module gw_iso_c
