#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module derivative_mod

  use derivative_mod_base, only:  &
      derivative_t,               &
      derivative_stag_t,          &
      subcell_integration,        &
      subcell_dss_fluxes,         &
      subcell_div_fluxes,         &
      subcell_Laplace_fluxes,     &
      allocate_subcell_integration_matrix, &
      derivinit,                  &
      deriv_print,                &
      gradient,                   &
      gradient_wk,                &
      vorticity,                  &
      divergence,                 &
      interpolate_gll2fvm_corners,&
      interpolate_gll2fvm_points, &
      remap_phys2gll,             &
      v2pinit,                    &
      gradient_sphere_wk_testcov, &
      gradient_sphere_wk_testcontra, &
      ugradv_sphere,              &
      vorticity_sphere,           &
      vorticity_sphere_diag,      &
      curl_sphere,                &
      curl_sphere_wk_testcov,     &
      vlaplace_sphere_wk,         &
      element_boundary_integral,  &
      edge_flux_u_cg,             &
      limiter_optim_iter_full,    &
      laplace_sphere_wk,          &
      divergence_sphere_wk,       &
      gradient_sphere,            &
      divergence_sphere

  implicit none

end module derivative_mod
