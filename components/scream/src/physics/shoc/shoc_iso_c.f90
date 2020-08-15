module shoc_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to shoc fortran.
!

contains
  subroutine append_precision(string, prefix)

    character(kind=c_char, len=128), intent(inout) :: string
    character(*), intent(in) :: prefix
    real(kind=c_real) :: s

    write (string, '(a,i1,a1)') prefix, sizeof(s), C_NULL_CHAR
  end subroutine append_precision

  subroutine shoc_init_c(nlev, gravit, rair, rh2o, cpair, &
                         zvir, latvap, latice, karman) bind(c)
    use shoc, only: shoc_init, npbl

    integer(kind=c_int), value, intent(in) :: nlev ! number of levels

    real(kind=c_real), value, intent(in)  :: gravit ! gravity
    real(kind=c_real), value, intent(in)  :: rair   ! dry air gas constant
    real(kind=c_real), value, intent(in)  :: rh2o   ! water vapor gas constant
    real(kind=c_real), value, intent(in)  :: cpair  ! specific heat of dry air
    real(kind=c_real), value, intent(in)  :: zvir   ! rh2o/rair - 1
    real(kind=c_real), value, intent(in)  :: latvap ! latent heat of vaporization
    real(kind=c_real), value, intent(in)  :: latice ! latent heat of fusion
    real(kind=c_real), value, intent(in)  :: karman ! Von Karman's constant

    real(kind=c_real) :: pref_mid(nlev) ! unused values

    pref_mid = 0
    call shoc_init(nlev, gravit, rair, rh2o, cpair, &
                   zvir, latvap, latice, karman, &
                   pref_mid, nlev, 1)
    npbl = nlev ! set pbl layer explicitly so we don't need pref_mid.
  end subroutine shoc_init_c

  subroutine shoc_main_c(shcol,nlev,nlevi,dtime,nadv,host_dx, host_dy, thv,  &
     zt_grid, zi_grid, pres, presi, pdel, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &
     wtracer_sfc, num_qtracers, w_field, exner,phis, host_dse, tke, thetal,  &
     qw, u_wind, v_wind, qtracers, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac, &
     pblh, shoc_mix, isotropy, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec,  &
     wqw_sec, wtke_sec, uw_sec, vw_sec, w3, wqls_sec, brunt, shoc_ql2) bind(C)
    use shoc, only : shoc_main

    integer(kind=c_int), value, intent(in) :: shcol, nlev, nlevi, num_qtracers, nadv
    real(kind=c_real), value, intent(in) :: dtime
    real(kind=c_real), intent(in), dimension(shcol) :: host_dx, host_dy
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: zt_grid
    real(kind=c_real), intent(in), dimension(shcol, nlevi) :: zi_grid
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: pres
    real(kind=c_real), intent(in), dimension(shcol, nlevi) :: presi
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: pdel, thv, w_field
    real(kind=c_real), intent(in), dimension(shcol) :: wthl_sfc, wqw_sfc, uw_sfc, vw_sfc
    real(kind=c_real), intent(in), dimension(shcol, num_qtracers) :: wtracer_sfc
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: exner
    real(kind=c_real), intent(in), dimension(shcol) :: phis

    real(kind=c_real), intent(inout), dimension(shcol, nlev) :: host_dse, tke, &
       thetal, qw, u_wind, v_wind, wthv_sec
    real(kind=c_real), intent(inout), dimension(shcol, nlev, num_qtracers) :: qtracers
    real(kind=c_real), intent(inout), dimension(shcol, nlev) :: tk, tkh

    real(kind=c_real), intent(out), dimension(shcol, nlev) :: shoc_cldfrac, shoc_ql
    real(kind=c_real), intent(out), dimension(shcol) :: pblh
    real(kind=c_real), intent(out), dimension(shcol, nlev) :: shoc_mix, w_sec
    real(kind=c_real), intent(out), dimension(shcol, nlevi) :: thl_sec, qw_sec, &
       qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3
    real(kind=c_real), intent(out), dimension(shcol, nlev) :: wqls_sec, isotropy, &
         brunt,shoc_ql2

    call shoc_main(shcol, nlev, nlevi, dtime, nadv, host_dx, host_dy, thv,   &
     zt_grid, zi_grid, pres, presi, pdel, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &
     wtracer_sfc, num_qtracers, w_field, exner, phis, host_dse, tke, thetal, &
     qw, u_wind, v_wind, qtracers, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac, &
     pblh, shoc_mix, isotropy, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec,  &
     wqw_sec, wtke_sec, uw_sec, vw_sec, w3, wqls_sec, brunt,shoc_ql2 )
  end subroutine shoc_main_c

  subroutine shoc_use_cxx_c(arg_use_cxx) bind(C)
    use shoc, only: use_cxx

    logical(kind=c_bool), value, intent(in) :: arg_use_cxx

    use_cxx = arg_use_cxx
  end subroutine shoc_use_cxx_c

  subroutine shoc_grid_c(shcol,nlev,nlevi,zt_grid,zi_grid,pdel,dz_zt,dz_zi,rho_zt) bind (C)
    use shoc, only: shoc_grid

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: zi_grid(shcol,nlevi)
    real(kind=c_real), intent(in) :: pdel(shcol,nlev)

    real(kind=c_real), intent(out) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(out) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(out) :: rho_zt(shcol,nlev)

    call shoc_grid(shcol,nlev,nlevi,zt_grid,zi_grid,pdel,dz_zt,dz_zi,rho_zt)

  end subroutine shoc_grid_c

  subroutine calc_shoc_varorcovar_c(&
       shcol,nlev,nlevi,tunefac,&                ! Input
       isotropy_zi,tkh_zi,dz_zi,invar1,invar2,&  ! Input
       varorcovar)bind (C)                       ! Input/Output

      use shoc, only: calc_shoc_varorcovar

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

      call calc_shoc_varorcovar(&
           shcol,nlev,nlevi,tunefac,&               ! Input
           isotropy_zi,tkh_zi,dz_zi,invar1,invar2,& ! Input
           varorcovar)

  end subroutine calc_shoc_varorcovar_c

  subroutine integ_column_stability_c(nlev, shcol, dz_zt, pres, brunt, brunt_int) bind (C)
    use shoc, only: integ_column_stability

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: pres(shcol,nlev)
    real(kind=c_real), intent(in) :: brunt(shcol,nlev)

    real(kind=c_real), intent(out) :: brunt_int(shcol)

    call integ_column_stability(nlev, shcol, dz_zt, pres, brunt, brunt_int)

  end subroutine integ_column_stability_c

  subroutine compute_shr_prod_c(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm) bind (C)
    use shoc, only: compute_shr_prod

    integer(kind=c_int), intent(in), value :: nlevi
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: u_wind(shcol,nlev)
    real(kind=c_real), intent(in) :: v_wind(shcol,nlev)

    real(kind=c_real), intent(out) :: sterm(shcol,nlevi)

    call compute_shr_prod(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm)

  end subroutine compute_shr_prod_c

  subroutine isotropic_ts_c(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy) bind (C)
    use shoc, only: isotropic_ts

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: brunt_int(shcol)
    real(kind=c_real), intent(in) :: tke(shcol,nlev)
    real(kind=c_real), intent(in) :: a_diss(shcol,nlev)
    real(kind=c_real), intent(in) :: brunt(shcol,nlev)

    real(kind=c_real), intent(out) :: isotropy(shcol,nlev)

    call isotropic_ts(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy)

  end subroutine isotropic_ts_c

  subroutine adv_sgs_tke_c(nlev, shcol, dtime, shoc_mix, wthv_sec, &
                           sterm_zt, tk, tke, a_diss) bind (C)
    use shoc, only: adv_sgs_tke

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in), value :: dtime
    real(kind=c_real), intent(in) :: shoc_mix(shcol,nlev)
    real(kind=c_real), intent(in) :: wthv_sec(shcol,nlev)
    real(kind=c_real), intent(in) :: sterm_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: tk(shcol,nlev)

    real(kind=c_real), intent(inout) :: tke(shcol,nlev)

    real(kind=c_real), intent(out) :: a_diss(shcol,nlev)

    call adv_sgs_tke(nlev, shcol, dtime, shoc_mix, wthv_sec, &
                     sterm_zt, tk, tke, a_diss)

  end subroutine adv_sgs_tke_c

  subroutine eddy_diffusivities_c(nlev, shcol, obklen, pblh, zt_grid, &
                          shoc_mix, sterm_zt, isotropy, tke, tkh, tk) bind (C)
    use shoc, only: eddy_diffusivities

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: obklen(shcol)
    real(kind=c_real), intent(in) :: pblh(shcol)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: shoc_mix(shcol,nlev)
    real(kind=c_real), intent(in) :: sterm_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: isotropy(shcol,nlev)
    real(kind=c_real), intent(in) :: tke(shcol,nlev)

    real(kind=c_real), intent(out) :: tkh(shcol,nlev)
    real(kind=c_real), intent(out) :: tk(shcol,nlev)

    call eddy_diffusivities(nlev, shcol, obklen, pblh, zt_grid, &
     shoc_mix, sterm_zt, isotropy, tke, tkh, tk)

  end subroutine eddy_diffusivities_c

  subroutine calc_shoc_vertflux_c(shcol, nlev, nlevi, tkh_zi, dz_zi, invar, vertflux) bind (C)
    use shoc, only: calc_shoc_vertflux

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in) :: tkh_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: invar(shcol,nlev)

    real(kind=c_real), intent(inout) :: vertflux(shcol,nlevi)

    call calc_shoc_vertflux(shcol, nlev, nlevi, tkh_zi, dz_zi, invar, vertflux)

  end subroutine calc_shoc_vertflux_c

  subroutine shoc_diag_second_moments_srf_c(shcol, wthl, uw, vw, ustar2, wstar) bind(C)
   use shoc, only: diag_second_moments_srf

   ! argmens
   integer(kind=c_int), value, intent(in) :: shcol
   real(kind=c_real), intent(in)  :: wthl(shcol), uw(shcol), vw(shcol)
   real(kind=c_real), intent(out) :: ustar2(shcol), wstar(shcol)

   call diag_second_moments_srf(shcol, wthl, uw, vw, ustar2, wstar)
 end subroutine shoc_diag_second_moments_srf_c

end module shoc_iso_c
