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

void gwd_compute_tendencies_from_stress_divergence_c(Int ncol, Int ngwv, bool do_taper, Real dt, Real effgw, Int* tend_level, Real* lat, Real* dpm, Real* rdpm, Real* c, Real* ubm, Real* t, Real* nm, Real* xv, Real* yv, Real* tau, Real* gwut, Real* utgw, Real* vtgw);

void gw_init_c(Int pver_in, Int pgwv_in, Real dc_in, Real* cref_in, bool do_molec_diff_in, bool tau_0_ubc_in, Int nbot_molec_in, Int ktop_in, Int kbotbg_in, Real fcrit2_in, Real kwv_in, Real gravit_in, Real rair_in, Real* alpha_in);

} // extern "C" : end _c decls

// Wrapper around gw_init
void gw_init(
  Int pver_in = 0,
  Int pgwv_in = 0,
  Real dc_in = 0.,
  Real* cref_in = nullptr,
  bool do_molec_diff_in = false,
  bool tau_0_ubc_in = false,
  Int nbot_molec_in = 0,
  Int ktop_in = 0,
  Int kbotbg_in = 0,
  Real fcrit2_in = 0.,
  Real kwv_in = 0.,
  Real gravit_in = 0.,
  Real rair_in = 0.,
  Real* alpha_in = nullptr)
{
  gw_init_c(pver_in, pgwv_in, dc_in, cref_in, do_molec_diff_in,tau_0_ubc_in, nbot_molec_in, ktop_in, kbotbg_in, fcrit2_in, kwv_in, gravit_in, rair_in, alpha_in);
}

void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d)
{
  gw_init(d.pver, d.pgwv);
  d.transpose<ekat::TransposeDirection::c2f>();
  gwd_compute_tendencies_from_stress_divergence_c(d.ncol, d.ngwv, d.do_taper, d.dt, d.effgw, d.tend_level, d.lat, d.dpm, d.rdpm, d.c, d.ubm, d.t, d.nm, d.xv, d.yv, d.tau, d.gwut, d.utgw, d.vtgw);
  d.transpose<ekat::TransposeDirection::f2c>();
}

// end _c impls

} // namespace gw
} // namespace scream
