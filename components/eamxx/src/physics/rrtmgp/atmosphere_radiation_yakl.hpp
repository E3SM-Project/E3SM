#ifndef ATMOSPHERE_RADIATION_YAKL_HPP
#define ATMOSPHERE_RADIATION_YAKL_HPP

#include "cpp/rrtmgp_const.h"

namespace scream {

void yakl_init ();

void compute_diffuse_flux(const int nswbands, const int nlay, const int ncol,
                          const real3d& sw_bnd_flux_dif,
                          const real3d& sw_bnd_flux_dn,
                          const real3d& sw_bnd_flux_dir);

void yakl_finalize();

} // namespace scream

#endif // ATMOSPHERE_RADIATION_YAKL_HPP
