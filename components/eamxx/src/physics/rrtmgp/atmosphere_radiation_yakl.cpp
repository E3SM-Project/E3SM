#include "atmosphere_radiation_yakl.hpp"
#include "YAKL.h"

namespace scream {

void yakl_init ()
{
  // Initialize yakl
  if(!yakl::isInitialized()) { yakl::init(); }
}

// Compute diffuse flux as difference between total and direct; use YAKL parallel_for here because these are YAKL objects
void compute_diffuse_flux(const int nswbands, const int nlay, const int ncol,
                          const real3d& sw_bnd_flux_dif,
                          const real3d& sw_bnd_flux_dn,
                          const real3d& sw_bnd_flux_dir)
{
  parallel_for(Bounds<3>(nswbands,nlay+1,ncol), YAKL_LAMBDA(int ibnd, int ilev, int icol) {
    sw_bnd_flux_dif(icol,ilev,ibnd) = sw_bnd_flux_dn(icol,ilev,ibnd) - sw_bnd_flux_dir(icol,ilev,ibnd);
  });
}

void yakl_finalize()
{
  // Finalize YAKL
  yakl::finalize();
}

}  // namespace scream
