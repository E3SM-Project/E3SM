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

  subroutine shoc_diag_second_moments_srf_f(shcol, wthl, uw, vw, ustar2, wstar) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: shcol

    ! arguments
    real(kind=c_real), intent(in) :: wthl(shcol)
    real(kind=c_real), intent(in) :: uw(shcol)
    real(kind=c_real), intent(in) :: vw(shcol)
    real(kind=c_real), intent(out) :: ustar2(shcol)
    real(kind=c_real), intent(out) :: wstar(shcol)

  end subroutine shoc_diag_second_moments_srf_f

  subroutine shoc_diag_second_moments_ubycond_f(shcol, thl, qw, wthl, wqw, qwthl, uw, vw, wtke) bind(C)
    use iso_c_binding

    ! argmens
    integer(kind=c_int), value, intent(in) :: shcol
    real(kind=c_real), intent(out)  :: thl(shcol), qw(shcol), qwthl(shcol),wthl(shcol),wqw(shcol), &
         uw(shcol), vw(shcol), wtke(shcol)

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

    real(kind=c_real) , intent(in), dimension(ncol, km1) :: x1, y1
    real(kind=c_real) , intent(in), dimension(ncol, km2) :: x2
    real(kind=c_real) , intent(out), dimension(ncol, km2) :: y2
    integer(kind=c_int) , value, intent(in) :: km1, km2, ncol
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

end interface

end module shoc_iso_f
