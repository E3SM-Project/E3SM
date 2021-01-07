module shoc_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from shoc fortran to scream c++.
!

interface

  subroutine calc_shoc_varorcovar_f(shcol, nlev, nlevi, tunefac, isotropy_zi, tkh_zi, dz_zi, invar1, invar2, varorcovar) bind (C)
    use iso_c_binding

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in), value :: tunefac
    real(kind=c_real), intent(in) :: isotropy_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: tkh_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: invar1(shcol,nlev)
    real(kind=c_real), intent(in) :: invar2(shcol,nlev)

    real(kind=c_real), intent(inout) :: varorcovar(shcol,nlevi)

  end subroutine calc_shoc_varorcovar_f

  subroutine calc_shoc_vertflux_f(shcol, nlev, nlevi, tkh_zi, dz_zi, invar, vertflux) bind (C)
    use iso_c_binding

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in) :: tkh_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: invar(shcol,nlev)

    real(kind=c_real), intent(inout) :: vertflux(shcol,nlevi)

  end subroutine calc_shoc_vertflux_f

  subroutine shoc_diag_second_moments_srf_f(shcol, wthl_sfc, uw_sfc, vw_sfc, ustar2, wstar) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: shcol

    ! arguments
    real(kind=c_real), intent(in) :: wthl_sfc(shcol)
    real(kind=c_real), intent(in) :: uw_sfc(shcol)
    real(kind=c_real), intent(in) :: vw_sfc(shcol)
    real(kind=c_real), intent(out) :: ustar2(shcol)
    real(kind=c_real), intent(out) :: wstar(shcol)

  end subroutine shoc_diag_second_moments_srf_f

  subroutine shoc_diag_second_moments_ubycond_f(shcol, thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec) bind(C)
    use iso_c_binding

    ! argmens
    integer(kind=c_int), value, intent(in) :: shcol
    real(kind=c_real), intent(out)  :: thl_sec(shcol), qw_sec(shcol), qwthl_sec(shcol),wthl_sec(shcol),wqw_sec(shcol), &
         uw_sec(shcol), vw_sec(shcol), wtke_sec(shcol)

  end subroutine shoc_diag_second_moments_ubycond_f

  subroutine update_host_dse_f(shcol, nlev, thlm, shoc_ql, exner, zt_grid, &
       phis, host_dse) bind (C)
    use iso_c_binding

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    real(kind=c_real), intent(in) :: thlm(shcol,nlev)
    real(kind=c_real), intent(in) :: shoc_ql(shcol,nlev)
    real(kind=c_real), intent(in) :: exner(shcol,nlev)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: phis(shcol)

    real(kind=c_real), intent(out) :: host_dse(shcol,nlev)

  end subroutine update_host_dse_f

  subroutine compute_diag_third_shoc_moment_f(shcol, nlev, nlevi, w_sec, thl_sec, &
                                              wthl_sec, tke, dz_zt, &
                                              dz_zi, isotropy_zi, &
                                              brunt_zi, w_sec_zi, thetal_zi, &
                                              w3) bind(C)
    use iso_c_binding

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in) :: w_sec(shcol,nlev)
    real(kind=c_real), intent(in) :: thl_sec(shcol,nlevi)
    real(kind=c_real), intent(in) :: wthl_sec(shcol,nlevi)
    real(kind=c_real), intent(in) :: tke(shcol,nlev)
    real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: isotropy_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: brunt_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: w_sec_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: thetal_zi(shcol,nlevi)

    real(kind=c_real), intent(out) :: w3(shcol,nlevi)

  end subroutine compute_diag_third_shoc_moment_f

  subroutine check_tke_f(shcol, nlev, tke) bind(C)

    use iso_c_binding

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev

    real(kind=c_real), intent(inout) :: tke(shcol,nlev)

  end subroutine check_tke_f

  subroutine shoc_pblintd_init_pot_f(shcol, nlev, thl, ql, q, thv) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: shcol, nlev
    real(kind=c_real), intent(in)  :: thl(shcol, nlev), ql(shcol, nlev), q(shcol,nlev)
    real(kind=c_real), intent(out) :: thv(shcol, nlev)

  end subroutine shoc_pblintd_init_pot_f

  subroutine compute_shoc_mix_shoc_length_f(nlev,shcol,tke,brunt,tscale,&
       zt_grid,l_inf,shoc_mix) bind (C)
    use iso_c_binding

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: tke(shcol,nlev)
    real(kind=c_real), intent(in) :: brunt(shcol,nlev)
    real(kind=c_real), intent(in) :: tscale(shcol)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: l_inf(shcol)

    real(kind=c_real), intent(out) :: shoc_mix(shcol,nlev)

  end subroutine compute_shoc_mix_shoc_length_f

  subroutine linear_interp_f(x1, x2, y1, y2, km1, km2, ncol, minthresh) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: km1, km2, ncol
    real(kind=c_real) , intent(in), dimension(ncol, km1) :: x1, y1
    real(kind=c_real) , intent(in), dimension(ncol, km2) :: x2
    real(kind=c_real) , intent(out), dimension(ncol, km2) :: y2
    real(kind=c_real) , value, intent(in) :: minthresh
  end subroutine linear_interp_f

subroutine clipping_diag_third_shoc_moments_f(nlevi,shcol,w_sec_zi,w3) bind (C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: nlevi
  integer(kind=c_int), intent(in), value :: shcol
  real(kind=c_real), intent(in) :: w_sec_zi(shcol,nlevi)

  real(kind=c_real), intent(inout) :: w3(shcol,nlevi)

end subroutine clipping_diag_third_shoc_moments_f

subroutine shoc_energy_integrals_f(shcol, nlev, host_dse, pdel,&
                                   rtm, rcm, u_wind, v_wind,&
                                   se_int, ke_int, wv_int, wl_int) bind (C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: shcol
  integer(kind=c_int), intent(in), value :: nlev
  real(kind=c_real), intent(in) :: host_dse(shcol,nlev)
  real(kind=c_real), intent(in) :: pdel(shcol,nlev)
  real(kind=c_real), intent(in) :: rtm(shcol,nlev)
  real(kind=c_real), intent(in) :: rcm(shcol,nlev)
  real(kind=c_real), intent(in) :: u_wind(shcol,nlev)
  real(kind=c_real), intent(in) :: v_wind(shcol,nlev)

  real(kind=c_real), intent(out) :: se_int(shcol)
  real(kind=c_real), intent(out) :: ke_int(shcol)
  real(kind=c_real), intent(out) :: wv_int(shcol)
  real(kind=c_real), intent(out) :: wl_int(shcol)

end subroutine shoc_energy_integrals_f

subroutine diag_second_moments_lbycond_f(shcol, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, ustar2, wstar, &
  wthl_sec, wqw_sec, uw_sec, vw_sec, wtke_sec, thl_sec, qw_sec, qwthl_sec) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: shcol
  real(kind=c_real) , intent(in), dimension(shcol) :: wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, ustar2, wstar
  real(kind=c_real) , intent(out), dimension(shcol) :: wthl_sec, wqw_sec, uw_sec, vw_sec, wtke_sec, thl_sec, qw_sec, qwthl_sec
end subroutine diag_second_moments_lbycond_f

subroutine diag_second_moments_f(shcol, nlev, nlevi, thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk, dz_zi, zt_grid, zi_grid, &
                                 shoc_mix, thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
  real(kind=c_real) , intent(in), dimension(shcol, nlev) :: thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk, zt_grid, shoc_mix
  real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: dz_zi, zi_grid
  real(kind=c_real) , intent(inout), dimension(shcol, nlevi) :: thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec
  real(kind=c_real) , intent(out), dimension(shcol, nlev) :: w_sec
end subroutine diag_second_moments_f

  subroutine diag_second_shoc_moments_f(shcol, nlev, nlevi, thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk, dz_zi, zt_grid, &
                                        zi_grid, shoc_mix, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, thl_sec, qw_sec, wthl_sec, wqw_sec, &
                                        qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: thetal, qw, u_wind, v_wind, tke, isotropy, tkh, tk, zt_grid, shoc_mix
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: dz_zi, zi_grid
    real(kind=c_real) , intent(in), dimension(shcol) :: wthl_sfc, wqw_sfc, uw_sfc, vw_sfc
    real(kind=c_real) , intent(out), dimension(shcol, nlevi) :: thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, uw_sec, vw_sec, wtke_sec
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: w_sec
  end subroutine diag_second_shoc_moments_f

subroutine compute_brunt_shoc_length_f(nlev, nlevi, shcol, dz_zt, thv, thv_zi, brunt) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: nlev, nlevi, shcol
  real(kind=c_real) , intent(in), dimension(shcol, nlev) :: dz_zt, thv
  real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: thv_zi
  real(kind=c_real) , intent(out), dimension(shcol, nlev) :: brunt

end subroutine compute_brunt_shoc_length_f

subroutine compute_l_inf_shoc_length_f(nlev,shcol,zt_grid,dz_zt,tke,l_inf) bind (C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: nlev
  integer(kind=c_int), intent(in), value :: shcol

  real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
  real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
  real(kind=c_real), intent(in) :: tke(shcol,nlev)

  real(kind=c_real), intent(out) :: l_inf(shcol)

end subroutine compute_l_inf_shoc_length_f

subroutine check_length_scale_shoc_length_f(nlev, shcol, host_dx, host_dy, shoc_mix) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: nlev, shcol
  real(kind=c_real) , intent(in), dimension(shcol) :: host_dx, host_dy
  real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: shoc_mix

end subroutine check_length_scale_shoc_length_f

subroutine compute_conv_vel_shoc_length_f(nlev,shcol,pblh,zt_grid,dz_zt,&
                                          thv,wthv_sec,conv_vel) bind (C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: nlev
  integer(kind=c_int), intent(in), value :: shcol
  real(kind=c_real), intent(in) :: pblh(shcol)
  real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
  real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
  real(kind=c_real), intent(in) :: thv(shcol,nlev)
  real(kind=c_real), intent(in) :: wthv_sec(shcol,nlev)

  real(kind=c_real), intent(out) :: conv_vel(shcol)

end subroutine compute_conv_vel_shoc_length_f

subroutine shoc_diag_obklen_f(shcol, uw_sfc, vw_sfc, wthl_sfc, wqw_sfc, thl_sfc, cldliq_sfc, qv_sfc, ustar, kbfs, obklen) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: shcol
  real(kind=c_real) , intent(in), dimension(shcol) :: uw_sfc, vw_sfc, wthl_sfc, wqw_sfc, thl_sfc, cldliq_sfc, qv_sfc
  real(kind=c_real) , intent(out), dimension(shcol) :: ustar, kbfs, obklen

end subroutine shoc_diag_obklen_f

subroutine shoc_pblintd_cldcheck_f(shcol, nlev, nlevi, zi, cldn, pblh) bind(C)
  use iso_c_binding

  integer(kind=c_int), value, intent(in) :: shcol, nlev, nlevi
  real(kind=c_real), intent(in) :: zi(shcol, nlevi), cldn(shcol, nlev)
  real(kind=c_real), intent(inout) :: pblh(shcol)
end subroutine shoc_pblintd_cldcheck_f

subroutine compute_conv_time_shoc_length_f(shcol,pblh,conv_vel,tscale) bind (C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: shcol
  real(kind=c_real), intent(in) :: pblh(shcol)
  real(kind=c_real), intent(inout) :: conv_vel(shcol)
  real(kind=c_real), intent(out) :: tscale(shcol)

end subroutine compute_conv_time_shoc_length_f

subroutine shoc_length_f(shcol, nlev, nlevi, host_dx, host_dy, pblh, &
                         tke, zt_grid, zi_grid, dz_zt, dz_zi, wthv_sec, thetal, thv, &
                         brunt, shoc_mix) bind (C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: shcol
  integer(kind=c_int), intent(in), value :: nlev
  integer(kind=c_int), intent(in), value :: nlevi

  real(kind=c_real), intent(in) :: host_dx(shcol)
  real(kind=c_real), intent(in) :: host_dy(shcol)
  real(kind=c_real), intent(in) :: pblh(shcol)

  real(kind=c_real), intent(in) :: tke(shcol,nlev)
  real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
  real(kind=c_real), intent(in) :: zi_grid(shcol,nlevi)
  real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
  real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
  real(kind=c_real), intent(in) :: wthv_sec(shcol,nlev)
  real(kind=c_real), intent(in) :: thetal(shcol,nlev)
  real(kind=c_real), intent(in) :: thv(shcol,nlev)

  real(kind=c_real), intent(out) :: brunt(shcol,nlev)
  real(kind=c_real), intent(out) :: shoc_mix(shcol,nlev)

end subroutine shoc_length_f

subroutine shoc_energy_fixer_f(shcol, nlev, nlevi, dtime, nadv, zt_grid, zi_grid,&
                               se_b, ke_b, wv_b, wl_b, se_a, ke_a, wv_a, wl_a,&
                               wthl_sfc, wqw_sfc, rho_zt, tke, pint, host_dse) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, nadv
  real(kind=c_real) , value, intent(in) :: dtime
  real(kind=c_real) , intent(in), dimension(shcol, nlev) :: zt_grid, rho_zt, tke
  real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: zi_grid, pint
  real(kind=c_real) , intent(in), dimension(shcol) :: se_b, ke_b, wv_b, wl_b, se_a, ke_a, wv_a, wl_a, wthl_sfc, wqw_sfc
  real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: host_dse

end subroutine shoc_energy_fixer_f

  subroutine compute_shoc_vapor_f(shcol, nlev, qw, ql, qv) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: qw, ql
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: qv
  end subroutine compute_shoc_vapor_f

  subroutine update_prognostics_implicit_f(shcol, nlev, nlevi, num_tracer, dtime, dz_zt, dz_zi, rho_zt, zt_grid, zi_grid, tk, tkh, uw_sfc, vw_sfc, wthl_sfc, wqw_sfc, wtracer_sfc, thetal, qw, tracer, tke, u_wind, v_wind) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, num_tracer
    real(kind=c_real) , value, intent(in) :: dtime
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: dz_zt, rho_zt, zt_grid, tk, tkh
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: dz_zi, zi_grid
    real(kind=c_real) , intent(in), dimension(shcol) :: uw_sfc, vw_sfc, wthl_sfc, wqw_sfc
    real(kind=c_real) , intent(in), dimension(shcol, num_tracer) :: wtracer_sfc
    real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: thetal, qw, tke, u_wind, v_wind
    real(kind=c_real) , intent(inout), dimension(shcol, nlev, num_tracer) :: tracer
  end subroutine update_prognostics_implicit_f

subroutine diag_third_shoc_moments_f(shcol, nlev, nlevi, w_sec, thl_sec, wthl_sec, isotropy, brunt,&
                                     thetal, tke, dz_zt, dz_zi, zt_grid, zi_grid, w3) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
  real(kind=c_real) , intent(in), dimension(shcol, nlev) :: w_sec, isotropy, brunt, thetal, tke, dz_zt, zt_grid
  real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: thl_sec, wthl_sec, dz_zi, zi_grid
  real(kind=c_real) , intent(out), dimension(shcol, nlevi) :: w3
end subroutine diag_third_shoc_moments_f

  subroutine adv_sgs_tke_f(nlev, shcol, dtime, shoc_mix, wthv_sec, sterm_zt, tk, tke, a_diss) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: nlev, shcol
    real(kind=c_real) , value, intent(in) :: dtime
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: shoc_mix, wthv_sec, sterm_zt, tk
    real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: tke
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: a_diss
  end subroutine adv_sgs_tke_f

subroutine shoc_assumed_pdf_f(shcol, nlev, nlevi, thetal, qw, w_field, thl_sec, qw_sec,&
                              wthl_sec, w_sec, wqw_sec, qwthl_sec, w3, pres, zt_grid,&
                              zi_grid, shoc_cldfrac, shoc_ql, wqls, wthv_sec, shoc_ql2) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
  real(kind=c_real) , intent(in), dimension(shcol, nlev) :: thetal, qw, w_field, w_sec, pres, zt_grid
  real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: thl_sec, qw_sec, wthl_sec, wqw_sec, qwthl_sec, w3, zi_grid
  real(kind=c_real) , intent(out), dimension(shcol, nlev) :: shoc_cldfrac, shoc_ql, wqls, wthv_sec, shoc_ql2
end subroutine shoc_assumed_pdf_f

subroutine compute_shr_prod_f(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm) bind(C)
  use iso_c_binding

  integer(kind=c_int) , value, intent(in) :: nlevi, nlev, shcol
  real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: dz_zi
  real(kind=c_real) , intent(in), dimension(shcol, nlev) :: u_wind, v_wind
  real(kind=c_real) , intent(out), dimension(shcol, nlevi) :: sterm
end subroutine compute_shr_prod_f

subroutine compute_tmpi_f(nlevi, shcol, dtime, rho_zi, dz_zi, tmpi) bind(C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: nlevi
  integer(kind=c_int), intent(in), value :: shcol
  real(kind=c_real), intent(in), value :: dtime
  real(kind=c_real), intent(in) :: rho_zi(shcol,nlevi)
  real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)

  real(kind=c_real), intent(out) :: tmpi(shcol,nlevi)
end subroutine compute_tmpi_f

subroutine integ_column_stability_f(nlev, shcol, dz_zt, pres, brunt, brunt_int) bind(C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: nlev
  integer(kind=c_int), intent(in), value :: shcol
  real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
  real(kind=c_real), intent(in) :: pres(shcol,nlev)
  real(kind=c_real), intent(in) :: brunt(shcol,nlev)
  real(kind=c_real), intent(out) :: brunt_int(shcol)

end subroutine integ_column_stability_f

subroutine dp_inverse_f(nlev, shcol, rho_zt, dz_zt, rdp_zt) bind(C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: nlev
  integer(kind=c_int), intent(in), value :: shcol
  real(kind=c_real), intent(in) :: rho_zt(shcol,nlev)
  real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
  real(kind=c_real), intent(out) :: rdp_zt(shcol,nlev)
end subroutine dp_inverse_f

  subroutine shoc_main_f(shcol, nlev, nlevi, dtime, nadv, host_dx, host_dy, thv, zt_grid, zi_grid, pres, presi, pdel, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, wtracer_sfc, num_qtracers, w_field, exner, phis, host_dse, tke, thetal, qw, u_wind, v_wind, qtracers, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac, pblh, shoc_mix, isotropy, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3, wqls_sec, brunt, shoc_ql2) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, nadv, num_qtracers
    real(kind=c_real) , value, intent(in) :: dtime
    real(kind=c_real) , intent(in), dimension(shcol) :: host_dx, host_dy, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, phis
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: thv, zt_grid, pres, pdel, w_field, exner
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: zi_grid, presi
    real(kind=c_real) , intent(in), dimension(shcol, num_qtracers) :: wtracer_sfc
    real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: host_dse, tke, thetal, qw, u_wind, v_wind, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac
    real(kind=c_real) , intent(inout), dimension(shcol, nlev, num_qtracers) :: qtracers
    real(kind=c_real) , intent(out), dimension(shcol) :: pblh
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: shoc_mix, isotropy, w_sec, wqls_sec, brunt, shoc_ql2
    real(kind=c_real) , intent(out), dimension(shcol, nlevi) :: thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3
  end subroutine shoc_main_f

  subroutine isotropic_ts_f(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: nlev, shcol
    real(kind=c_real) , intent(in), dimension(shcol) :: brunt_int
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: tke, a_diss, brunt
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: isotropy
  end subroutine isotropic_ts_f

  subroutine pblintd_height_f(shcol, nlev, z, u, v, ustar, thv, thv_ref, pblh, rino, check) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: z, u, v, thv
    real(kind=c_real) , intent(in), dimension(shcol) :: ustar, thv_ref
    real(kind=c_real) , intent(out), dimension(shcol) :: pblh
    real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: rino
    logical(kind=c_bool) , intent(inout), dimension(shcol) :: check
  end subroutine pblintd_height_f

  subroutine vd_shoc_decomp_and_solve_f(shcol, nlev, nlevi, num_rhs, kv_term, tmpi, rdp_zt, dtime, flux, var) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, num_rhs
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: kv_term, tmpi
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: rdp_zt
    real(kind=c_real) , value, intent(in) :: dtime
    real(kind=c_real) , intent(in), dimension(shcol) :: flux
    real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: var
  end subroutine vd_shoc_decomp_and_solve_f

  subroutine pblintd_surf_temp_f(shcol, nlev, nlevi, z, ustar, obklen, kbfs, thv, tlv, pblh, check, rino) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: z, thv
    real(kind=c_real) , intent(in), dimension(shcol) :: ustar, obklen, kbfs
    real(kind=c_real) , intent(out), dimension(shcol) :: tlv
    real(kind=c_real) , intent(inout), dimension(shcol) :: pblh
    logical(kind=c_bool) , intent(inout), dimension(shcol) :: check
    real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: rino
  end subroutine pblintd_surf_temp_f

  subroutine pblintd_check_pblh_f(shcol, nlev, nlevi, npbl, z, ustar, check, pblh) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi, npbl
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: z
    real(kind=c_real) , intent(in), dimension(shcol) :: ustar
    logical(kind=c_bool) , intent(in), dimension(shcol) :: check
    real(kind=c_real) , intent(inout), dimension(shcol) :: pblh
  end subroutine pblintd_check_pblh_f

  subroutine pblintd_f(shcol, nlev, nlevi, z, zi, thl, ql, q, u, v, ustar, obklen, kbfs, cldn, pblh) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: z, thl, ql, q, u, v, cldn
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: zi
    real(kind=c_real) , intent(in), dimension(shcol) :: ustar, obklen, kbfs
    real(kind=c_real) , intent(out), dimension(shcol) :: pblh
  end subroutine pblintd_f

  subroutine shoc_grid_f(shcol, nlev, nlevi, zt_grid, zi_grid, pdel, dz_zt, dz_zi, rho_zt) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: zt_grid, pdel
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: zi_grid
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: dz_zt, rho_zt
    real(kind=c_real) , intent(out), dimension(shcol, nlevi) :: dz_zi
  end subroutine shoc_grid_f

  subroutine eddy_diffusivities_f(nlev, shcol, obklen, pblh, zt_grid, shoc_mix, sterm_zt, isotropy, tke, tkh, tk) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: nlev, shcol
    real(kind=c_real) , intent(in), dimension(shcol) :: obklen, pblh
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: zt_grid, shoc_mix, sterm_zt, isotropy, tke
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: tkh, tk
  end subroutine eddy_diffusivities_f

  subroutine shoc_tke_f(shcol, nlev, nlevi, dtime, wthv_sec, shoc_mix, dz_zi, dz_zt, pres, u_wind, v_wind, brunt, obklen, zt_grid, zi_grid, pblh, tke, tk, tkh, isotropy) bind(C)
    use iso_c_binding

    integer(kind=c_int) , value, intent(in) :: shcol, nlev, nlevi
    real(kind=c_real) , value, intent(in) :: dtime
    real(kind=c_real) , intent(in), dimension(shcol, nlev) :: wthv_sec, shoc_mix, dz_zt, pres, u_wind, v_wind, brunt, zt_grid
    real(kind=c_real) , intent(in), dimension(shcol, nlevi) :: dz_zi, zi_grid
    real(kind=c_real) , intent(in), dimension(shcol) :: obklen, pblh
    real(kind=c_real) , intent(inout), dimension(shcol, nlev) :: tke, tk, tkh
    real(kind=c_real) , intent(out), dimension(shcol, nlev) :: isotropy
  end subroutine shoc_tke_f
end interface

end module shoc_iso_f
