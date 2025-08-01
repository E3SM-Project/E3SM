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

  subroutine gwd_compute_tendencies_from_stress_divergence_c(ncol, do_taper, dt, effgw, tend_level, lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw) bind(C)
    use gw_common, only : gwd_compute_tendencies_from_stress_divergence, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: ncol
    logical(kind=c_bool) , value, intent(in) :: do_taper
    real(kind=c_real) , value, intent(in) :: dt, effgw
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol) :: lat, xv, yv
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: dpm, rdpm, ubm, t, nm
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(inout), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver, -pgwv:pgwv) :: gwut
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: utgw, vtgw

    call gwd_compute_tendencies_from_stress_divergence(ncol, pgwv, do_taper, dt, effgw, tend_level, lat, dpm, rdpm, c, ubm, t, nm, xv, yv, tau, gwut, utgw, vtgw)
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

  subroutine gwd_compute_stress_profiles_and_diffusivities_c(ncol, src_level, ubi, c, rhoi, ni, kvtt, t, ti, piln, tau) bind(C)
    use gw_common, only : gwd_compute_stress_profiles_and_diffusivities, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: ncol
    integer(kind=c_int) , intent(in), dimension(ncol) :: src_level
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: ubi, rhoi, ni, kvtt, ti, piln
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: t
    real(kind=c_real) , intent(inout), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau

    call gwd_compute_stress_profiles_and_diffusivities(ncol, pgwv, src_level, ubi, c, rhoi, ni, kvtt, t, ti, piln, tau)
  end subroutine gwd_compute_stress_profiles_and_diffusivities_c

  subroutine gwd_project_tau_c(ncol, tend_level, tau, ubi, c, xv, yv, taucd) bind(C)
    use gw_common, only : gwd_project_tau, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: ncol
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: ubi
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(in), dimension(ncol) :: xv, yv
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver, 4) :: taucd

    call gwd_project_tau(ncol, pgwv, tend_level, tau, ubi, c, xv, yv, taucd)
  end subroutine gwd_project_tau_c

  subroutine gwd_precalc_rhoi_c(pcnst, ncol, dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm, c, q, dse, egwdffi, qtgw, dttdf, dttke, ttgw) bind(C)
    use gw_common, only : gwd_precalc_rhoi, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: pcnst, ncol
    real(kind=c_real) , value, intent(in) :: dt
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: pmid, t, ubm, nm, rdpm, dse
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: pint
    real(kind=c_real) , intent(in), dimension(ncol, pver, -pgwv:pgwv) :: gwut
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(in), dimension(ncol, pver, pcnst) :: q
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: egwdffi
    real(kind=c_real) , intent(out), dimension(ncol, pver, pcnst) :: qtgw
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: dttdf, dttke, ttgw

    call gwd_precalc_rhoi(ncol, pgwv, dt, tend_level, pmid, pint, t, gwut, ubm, nm, rdpm, c, q, dse, egwdffi, qtgw, dttdf, dttke, ttgw)
  end subroutine gwd_precalc_rhoi_c

  subroutine gw_drag_prof_c(pcnst, ncol, src_level, tend_level, do_taper, dt, lat, t, ti, pmid, pint, dpm, rdpm, piln, rhoi, nm, ni, ubm, ubi, xv, yv, effgw, c, kvtt, q, dse, tau, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, gwut, dttdf, dttke) bind(C)
    use gw_common, only : gw_drag_prof, pver, pgwv

    integer(kind=c_int) , value, intent(in) :: pcnst, ncol
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
    real(kind=c_real) , intent(out), dimension(ncol, pver, -pgwv:pgwv) :: gwut

    call gw_drag_prof(ncol, pgwv, src_level, tend_level, do_taper, dt, lat, t, ti, pmid, pint, dpm, rdpm, piln, rhoi, nm, ni, ubm, ubi, xv, yv, effgw, c, kvtt, q, dse, tau, utgw, vtgw, ttgw, qtgw, taucd, egwdffi, gwut, dttdf, dttke)
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

  subroutine gw_front_gw_sources_c(ncol, kbot, frontgf, tau) bind(C)
    use gw_common, only : pver, pgwv
    use gw_front, only : gw_front_gw_sources

    integer(kind=c_int) , value, intent(in) :: ncol, kbot
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: frontgf
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau

    call gw_front_gw_sources(ncol, pgwv, kbot, frontgf, tau)
  end subroutine gw_front_gw_sources_c

  subroutine gw_cm_src_c(ncol, kbot, u, v, frontgf, src_level, tend_level, tau, ubm, ubi, xv, yv, c) bind(C)
    use gw_common, only : pver, pgwv
    use gw_front, only : gw_cm_src

    integer(kind=c_int) , value, intent(in) :: ncol, kbot
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: u, v
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: frontgf
    integer(kind=c_int) , intent(out), dimension(ncol) :: src_level, tend_level
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: ubm
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: ubi
    real(kind=c_real) , intent(out), dimension(ncol) :: xv, yv
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv) :: c

    call gw_cm_src(ncol, pgwv, kbot, u, v, frontgf, src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  end subroutine gw_cm_src_c

  subroutine gw_convect_init_c(maxh, maxuh, plev_src_wind, mfcc_in) bind(C)
    use gw_common, only : pgwv
    use gw_convect, only : gw_convect_init

    integer(kind=c_int) , value, intent(in) :: maxh, maxuh
    real(kind=c_real) , value, intent(in) :: plev_src_wind
    real(kind=c_real) , intent(in), dimension(maxh, -maxuh:maxuh, -pgwv:pgwv) :: mfcc_in

    character(len=128) :: errstring

    call gw_convect_init(plev_src_wind, mfcc_in, errstring)
  end subroutine gw_convect_init_c

  subroutine gw_convect_project_winds_c(ncol, u, v, xv, yv, ubm, ubi) bind(C)
    use gw_common, only : pver
    use gw_convect, only : gw_convect_project_winds

    integer(kind=c_int) , value, intent(in) :: ncol
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: u, v
    real(kind=c_real) , intent(out), dimension(ncol) :: xv, yv
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: ubm
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: ubi

    call gw_convect_project_winds(ncol, u, v, xv, yv, ubm, ubi)
  end subroutine gw_convect_project_winds_c

  subroutine gw_heating_depth_c(ncol, maxq0_conversion_factor, hdepth_scaling_factor, use_gw_convect_old, zm, netdt, mini, maxi, hdepth, maxq0_out, maxq0) bind(C)
    use gw_common, only : pver
    use gw_convect, only : gw_heating_depth

    integer(kind=c_int) , value, intent(in) :: ncol
    real(kind=c_real) , value, intent(in) :: maxq0_conversion_factor, hdepth_scaling_factor
    logical(kind=c_bool) , value, intent(in) :: use_gw_convect_old
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: zm
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: netdt
    integer(kind=c_int) , intent(out), dimension(ncol) :: mini, maxi
    real(kind=c_real) , intent(out), dimension(ncol) :: hdepth, maxq0_out, maxq0

    call gw_heating_depth(ncol, maxq0_conversion_factor, hdepth_scaling_factor, use_gw_convect_old, zm, netdt, mini, maxi, hdepth, maxq0_out, maxq0)
  end subroutine gw_heating_depth_c

  subroutine gw_storm_speed_c(ncol, storm_speed_min, ubm, mini, maxi, storm_speed, uh, umin, umax) bind(C)
    use gw_common, only : pver
    use gw_convect, only : gw_storm_speed

    integer(kind=c_int) , value, intent(in) :: ncol
    real(kind=c_real) , value, intent(in) :: storm_speed_min
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: ubm
    integer(kind=c_int) , intent(in), dimension(ncol) :: mini, maxi
    integer(kind=c_int) , intent(out), dimension(ncol) :: storm_speed
    real(kind=c_real) , intent(out), dimension(ncol) :: uh, umin, umax

    call gw_storm_speed(ncol, storm_speed_min, ubm, mini, maxi, storm_speed, uh, umin, umax)
  end subroutine gw_storm_speed_c

  subroutine gw_convect_gw_sources_c(ncol, lat, hdepth_min, hdepth, mini, maxi, netdt, uh, storm_speed, maxq0, umin, umax, tau) bind(C)
    use gw_common, only : pver, pgwv
    use gw_convect, only : gw_convect_gw_sources

    integer(kind=c_int) , value, intent(in) :: ncol
    real(kind=c_real) , intent(in), dimension(ncol) :: lat, hdepth, uh, maxq0, umin, umax
    real(kind=c_real) , value, intent(in) :: hdepth_min
    integer(kind=c_int) , intent(in), dimension(ncol) :: mini, maxi, storm_speed
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: netdt
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau

    call gw_convect_gw_sources(ncol, pgwv, lat, hdepth_min, hdepth, mini, maxi, netdt, uh, storm_speed, maxq0, umin, umax, tau)
  end subroutine gw_convect_gw_sources_c

  subroutine gw_beres_src_c(ncol, lat, u, v, netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth, maxq0_out, maxq0_conversion_factor, hdepth_scaling_factor, hdepth_min, storm_speed_min, use_gw_convect_old) bind(C)
    use gw_common, only : pver, pgwv
    use gw_convect, only : gw_beres_src

    integer(kind=c_int) , value, intent(in) :: ncol
    real(kind=c_real) , intent(in), dimension(ncol) :: lat
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: u, v, zm
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: netdt
    integer(kind=c_int) , intent(out), dimension(ncol) :: src_level, tend_level
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: ubm
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: ubi
    real(kind=c_real) , intent(out), dimension(ncol) :: xv, yv, hdepth, maxq0_out
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , value, intent(in) :: maxq0_conversion_factor, hdepth_scaling_factor, hdepth_min, storm_speed_min
    logical(kind=c_bool) , value, intent(in) :: use_gw_convect_old

    call gw_beres_src(ncol, pgwv, lat, u, v, netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth, maxq0_out, maxq0_conversion_factor, hdepth_scaling_factor, hdepth_min, storm_speed_min, use_gw_convect_old)
  end subroutine gw_beres_src_c

  subroutine gw_ediff_c(ncol, kbot, ktop, tend_level, gwut, ubm, nm, rho, dt, gravit, pmid, rdpm, c, egwdffi, decomp_ca, decomp_cc, decomp_dnom, decomp_ze) bind(C)
    use gw_common, only : pver, pgwv
    use gw_diffusion, only : gw_ediff
    use vdiff_lu_solver, only: lu_decomp

    integer(kind=c_int) , value, intent(in) :: ncol, kbot, ktop
    integer(kind=c_int) , intent(in), dimension(ncol) :: tend_level
    real(kind=c_real) , intent(in), dimension(ncol, pver, -pgwv:pgwv) :: gwut
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: ubm, nm, pmid, rdpm
    real(kind=c_real) , intent(in), dimension(ncol, pver+1) :: rho
    real(kind=c_real) , value, intent(in) :: dt, gravit
    real(kind=c_real) , intent(in), dimension(ncol, -pgwv:pgwv) :: c
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: egwdffi
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: decomp_ca, decomp_cc, decomp_dnom, decomp_ze

    type(lu_decomp) :: decomp

    call gw_ediff(ncol, pver, pgwv, kbot, ktop, tend_level, gwut, ubm, nm, rho, dt, gravit, pmid, rdpm, c, egwdffi, decomp)
    decomp_ca = decomp%ca
    decomp_cc = decomp%cc
    decomp_dnom = decomp%dnom
    decomp_ze = decomp%ze

    call decomp%finalize()

  end subroutine gw_ediff_c

  subroutine gw_diff_tend_c(ncol, kbot, ktop, q, dt, decomp_ca, decomp_cc, decomp_dnom, decomp_ze, dq) bind(C)
    use gw_common, only : pver
    use gw_diffusion, only : gw_diff_tend
    use vdiff_lu_solver, only: lu_decomp

    integer(kind=c_int) , value, intent(in) :: ncol, kbot, ktop
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: q
    real(kind=c_real) , value, intent(in) :: dt
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: decomp_ca, decomp_cc, decomp_dnom, decomp_ze
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: dq

    type(lu_decomp) :: decomp
    decomp = lu_decomp(ncol, pver)

    call gw_diff_tend(ncol, pver, kbot, ktop, q, dt, decomp, dq)

    call decomp%finalize()
  end subroutine gw_diff_tend_c

  subroutine gw_oro_init_c() bind(C)
    use gw_oro, only : gw_oro_init

    character(len=128) :: errstring

    call gw_oro_init(errstring)
  end subroutine gw_oro_init_c

  subroutine gw_oro_src_c(ncol, u, v, t, sgh, pmid, pint, dpm, zm, nm, src_level, tend_level, tau, ubm, ubi, xv, yv, c) bind(C)
    use gw_common, only : pver, pgwv
    use gw_oro, only : gw_oro_src

    integer(kind=c_int) , value, intent(in) :: ncol
    real(kind=c_real) , intent(in), dimension(ncol, pver) :: u, v, t, pmid, dpm, zm, nm
    real(kind=c_real) , intent(in), dimension(ncol) :: sgh
    real(kind=c_real) , intent(in), dimension(ncol, 0:pver) :: pint
    integer(kind=c_int) , intent(out), dimension(ncol) :: src_level, tend_level
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv, 0:pver) :: tau
    real(kind=c_real) , intent(out), dimension(ncol, pver) :: ubm
    real(kind=c_real) , intent(out), dimension(ncol, 0:pver) :: ubi
    real(kind=c_real) , intent(out), dimension(ncol) :: xv, yv
    real(kind=c_real) , intent(out), dimension(ncol, -pgwv:pgwv) :: c

    call gw_oro_src(ncol, u, v, t, sgh, pmid, pint, dpm, zm, nm, src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  end subroutine gw_oro_src_c
end module gw_iso_c
