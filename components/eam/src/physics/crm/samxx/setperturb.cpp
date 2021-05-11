#include "setperturb.h"
void setperturb() {
  // Add random noise near the surface to help turbulence develop
  // The random generator is seeded based on the global column id, 
  // which avoids a problematic sensitivity to pcols.
  int  constexpr perturb_num_layers  = 5;    // Number of levels to perturb
  real constexpr perturb_t_magnitude = 1.0;  // perturbation LSE amplitube [K]
  real factor_xy = 1. / (nx*ny);
  auto t_host     = t.createHostCopy();
  auto t0_host    = t0.createHostCopy();
  auto gcolp_host = gcolp.createHostCopy();
  // Apply random liquid static energy (LSE) perturbations
  for (int icrm = 0; icrm < ncrms; icrm++) {
    srand(gcolp_host(icrm)*1000);
    for (int k = 0; k < perturb_num_layers; k++) {
      // set perturb_k_scaling so that perturbation magnitude decreases with altitude
      real perturb_k_scaling = ((real)perturb_num_layers-k) / (real)perturb_num_layers;
      real t02 = 0;
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          // Generate a uniform random number in interval (0,1)
          real rand_perturb = static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );
          // onvert perturbation range from (0,1) to (-1,1)
          rand_perturb = 1.-2.*rand_perturb;
          // apply perturbation 
          t_host(k,j+offy_s,i+offx_s,icrm) = t_host(k,j+offy_s,i+offx_s,icrm) + rand_perturb * perturb_t_magnitude * perturb_k_scaling;
          // Calculate new average LSE for energy conservation scaling below
          t02 += t_host(k,j+offy_s,i+offx_s,icrm)*factor_xy;
        }
      }
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          t_host(k,j+offy_s,i+offx_s,icrm) = t_host(k,j+offy_s,i+offx_s,icrm) * t0_host(k,icrm) / t02;
        }
      }
    }
  }
  t_host.deep_copy_to(t);
}
