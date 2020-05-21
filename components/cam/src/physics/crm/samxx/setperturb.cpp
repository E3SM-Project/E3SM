
#include "setperturb.h"

void setperturb() {
  // Add random noise near the surface to help turbulence develop
  // This surboutine has been updated for SPCAM5 (Minghuai.Wang@pnnl.gov, April, 2012).
  // Now the random generator is seeded based on the global column id, which gets rid
  // of the dependence of the SPCAM results on pcols.
  // Walter Hannah - LLNL - Mar 2018
  int constexpr perturb_seed_scale  = 1000; // scaling value for setperturb() seed value (seed = gcol * perturb_seed_scale)
  int constexpr perturb_num_layers  = 5;    // Number of levels to perturb
  int constexpr perturb_t_magnitude = 1.0;  // perturbation LSE amplitube [K]

  real2d t02 = real2d("t02",perturb_num_layers,ncrms);

  real factor_xy = 1./(nx*ny);
  
  // Apply random liquid static energy (LSE) perturbations
  // do k = 1,perturb_num_layers
  //   do j = 1,ny
  //     do i = 1,nx
  //       do icrm = 1 , ncrms
  parallel_for( Bounds<4>(pergurb_num_layers,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    // set perturb_k_scaling so that perturbation magnitude decreases with altitude
    real perturb_k_scaling = ( perturb_num_layers+1.-k) / perturb_num_layers;
    // Get the random number
    yakl::Random rand(gcolp(icrm) * perturb_seed_scale + k*ny*nx*ncrms + j*nx*ncrms + i*ncrms + icrm);
    real rand_perturb = rand.genFP<real>( -1._fp , 1._fp );
    // apply perturbation 
    t(k,j,i,icrm) = t(k,j,i,icrm) + rand_perturb * perturb_t_magnitude * perturb_k_scaling;
    // Calculate new average LSE for energy conservation scaling below
    yakl::atomicAdd( t02(k,icrm) , t(k,j,i,icrm)*factor_xy );
  });

  // enforce energy conservation
  // do k = 1,perturb_num_layers
  //   do j = 1,ny
  //     do i = 1,nx
  //       do icrm = 1 , ncrms
  parallel_for( Bounds<4>(pergurb_num_layers,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    t(k,j,i,icrm) = t(k,j,i,icrm) * t0(k,icrm) / t02(k,icrm);
  });
}

