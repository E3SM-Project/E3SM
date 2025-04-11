#include "gw_test_data.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to gw fortran calls and vice versa
//

namespace scream {
namespace gw {

extern "C" {
void gwd_compute_tendencies_from_stress_divergence_c(Int ncol, Int pver, Int pgwv, Int ngwv, bool do_taper, Real dt, Real effgw, Int* tend_level, Real* lat, Real* dpm, Real* rdpm, Real* c, Real* ubm, Real* t, Real* nm, Real* xv, Real* yv, Real* tau, Real* gwut, Real* utgw, Real* vtgw);
} // extern "C" : end _c decls

void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d)
{
  gw_init();
  d.transpose<ekat::TransposeDirection::c2f>();
  gwd_compute_tendencies_from_stress_divergence_c(d.ncol, d.pver, d.pgwv, d.ngwv, d.do_taper, d.dt, d.effgw, d.tend_level, d.lat, d.dpm, d.rdpm, d.c, d.ubm, d.t, d.nm, d.xv, d.yv, d.tau, d.gwut, d.utgw, d.vtgw);
  d.transpose<ekat::TransposeDirection::f2c>();
}

// end _c impls

} // namespace gw
} // namespace scream
