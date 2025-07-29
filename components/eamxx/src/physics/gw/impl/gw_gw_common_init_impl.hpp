#ifndef GW_GW_COMMON_INIT_IMPL_HPP
#define GW_GW_COMMON_INIT_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_common_init. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::gw_common_init(
// Inputs
  const Int& pver_in,
  const Int& pgwv_in,
  const Real& dc_in,
  const uview_1d<const Real>& cref_in,
  const bool& orographic_only_in,
  const bool& do_molec_diff_in,
  const bool& tau_0_ubc_in,
  const Int& nbot_molec_in,
  const Int& ktop_in,
  const Int& kbotbg_in,
  const Real& fcrit2_in,
  const Real& kwv_in,
  const uview_1d<const Real>& alpha_in)
{
  s_common_init.initialized = true;
  s_common_init.pver = pver_in;
  s_common_init.pgwv = pgwv_in;
  s_common_init.dc = dc_in;
  s_common_init.cref = view_1d<Real>("cref", cref_in.size());
  Kokkos::deep_copy(s_common_init.cref, cref_in);
  s_common_init.orographic_only = orographic_only_in;
  s_common_init.do_molec_diff = do_molec_diff_in;
  s_common_init.tau_0_ubc = tau_0_ubc_in;
  s_common_init.nbot_molec = nbot_molec_in;
  s_common_init.ktop = ktop_in;
  s_common_init.kbotbg = kbotbg_in;
  s_common_init.fcrit2 = fcrit2_in;
  s_common_init.kwv = kwv_in;
  s_common_init.alpha = view_1d<Real>("alpha", alpha_in.size());
  Kokkos::deep_copy(s_common_init.alpha, alpha_in);
  s_common_init.effkwv = kwv_in * fcrit2_in;
  s_common_init.tndmax = orographic_only_in ? 500. / 86400. : 400. / 86400.;
}

} // namespace gw
} // namespace scream

#endif
