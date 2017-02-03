#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module derivative_mod
  use derivative_mod_base, only: &
    allocate_subcell_integration_matrix,&
    curl_sphere, &
    curl_sphere_wk_testcov,&
    derivative_t,&
    derivinit,&
    divergence,&
    divergence_sphere,&
    divergence_sphere_wk,&
    edge_flux_u_cg,&
    element_boundary_integral,&
    gradient,&
    gradient_sphere,&
    gradient_sphere_wk_testcontra,&
    gradient_sphere_wk_testcov,&
    gradient_wk,&
    laplace_sphere_wk,&
    limiter_optim_iter_full,&
    subcell_Laplace_fluxes,&
    subcell_div_fluxes,&
    subcell_dss_fluxes,&
    subcell_integration,&
    ugradv_sphere, vorticity_sphere,&
    vlaplace_sphere_wk,&
    vorticity,&
    vorticity_sphere_diag

  implicit none
end module derivative_mod
