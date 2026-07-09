#ifndef GW_GW_COMMON_INIT_IMPL_HPP
#define GW_GW_COMMON_INIT_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "share/util/eamxx_simple_linear_interp.hpp"

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
  static bool s_common_init_constructed = false;
  if (!s_common_init_constructed) {
    new (&s_common_init) GwCommonInit();
    s_common_init_constructed = true;
  }

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
  s_common_init.oroko2 = GWC::half * kwv_in;
  s_common_init.alpha = view_1d<Real>("alpha", alpha_in.size());
  Kokkos::deep_copy(s_common_init.alpha, alpha_in);
  s_common_init.effkwv = kwv_in * fcrit2_in;
  s_common_init.tndmax = orographic_only_in ? 500/GWC::sec_per_day : 400/GWC::sec_per_day;
}

/*
The version above was written for unit tests. The version below
is better suited for running the GW schemes in the full model
*/

template<typename S, typename D>
void Functions<S,D>::gw_common_init(
  // Inputs
  ekat::ParameterList& params,
  const Int& pver_in,
  const uview_1d<const Real>& pref_int,
  const bool& do_molec_diff_in,
  const Int& nbot_molec_in,
  const Int& ktop_in,
  const Real& kwv_in)
{
  // Intel compiler may not call constructors for inline static members of
  // template classes, leaving Kokkos view trackers zero-initialized (which
  // causes crashes on assignment). Use placement new to force-construct on
  // first call. Function-local statics are guaranteed to initialize correctly.
  static bool s_common_init_constructed = false;
  if (!s_common_init_constructed) {
    new (&s_common_init) GwCommonInit();
    s_common_init_constructed = true;
  }

  EKAT_REQUIRE_MSG(static_cast<int>(pref_int.size())==pver_in+1, "Error! pref_int size is incorrect.\n");

  using exe_space_t = typename KT::ExeSpace;

  s_common_init.initialized = true;
  s_common_init.pver = pver_in;

  s_common_init.use_gw_convect    = params.get<bool>("use_gw_convect", s_common_init.use_gw_convect);
  s_common_init.use_gw_frontal    = params.get<bool>("use_gw_frontal", s_common_init.use_gw_frontal);
  s_common_init.use_gw_orographic = params.get<bool>("use_gw_orographic", s_common_init.use_gw_orographic);
  s_common_init.pgwv              = params.get<int>("pgwv", s_common_init.pgwv);
  s_common_init.dc                = params.get<Real>("gw_dc", s_common_init.dc);
  s_common_init.tau_0_ubc         = params.get<bool>("tau_0_ubc", s_common_init.tau_0_ubc);
  s_common_init.fcrit2            = params.get<Real>("fcrit2", s_common_init.fcrit2);
  s_common_init.gw_orographic_eff = params.get<Real>("gw_orographic_eff", s_common_init.gw_orographic_eff);

  // Set phase speeds
  const int num_pgwv = s_common_init.pgwv * 2 + 1;
  const int pgwv = s_common_init.pgwv;
  const Real dc = s_common_init.dc;
  s_common_init.cref = view_1d<Real>("cref", num_pgwv);
  auto cref = s_common_init.cref;
  Kokkos::parallel_for("gw_common_init_cref",
    Kokkos::RangePolicy<exe_space_t>(0, num_pgwv),
    KOKKOS_LAMBDA(const int i) {
      cref(i) = dc * (i - pgwv);
    });

  s_common_init.orographic_only = !s_common_init.use_gw_convect && !s_common_init.use_gw_frontal;

  // calculate alpha (Newtonian cooling coefficients)
  s_common_init.alpha = view_1d<Real>("alpha", pver_in+1);
  if (s_common_init.orographic_only) {
    Kokkos::deep_copy(s_common_init.alpha,1e-6);
  } else {
    // pre-calculated newtonian damping:
    //     * convert alpha0 from 1/day to 1/s
    //     * ensure alpha0 is not smaller than 1e-6
    //     * convert alpha_pressure_mb to pa
    view_1d<Real> alpha0_per_sec("alpha0_per_sec", GWC::nalph);
    view_1d<Real> alpha_pressure_pa("alpha_pressure_pa", GWC::nalph);
    auto alpha0_per_sec_h   = Kokkos::create_mirror_view(alpha0_per_sec);
    auto alpha_pressure_pa_h = Kokkos::create_mirror_view(alpha_pressure_pa);
    for (int k = 0; k < GWC::nalph; ++k) {
      alpha0_per_sec_h(k)    = std::max(GWC::alpha0[k] / 86400.0, 1.0e-6);
      alpha_pressure_pa_h(k) = GWC::alpha_pressure_mb[k] * 1.0e2;
    }
    Kokkos::deep_copy(alpha0_per_sec,   alpha0_per_sec_h);
    Kokkos::deep_copy(alpha_pressure_pa, alpha_pressure_pa_h);
    simple_linear_interp(alpha_pressure_pa, alpha0_per_sec, pref_int, s_common_init.alpha);
  }

  // set bottom index for background spectrum (small serial scan; do on host)
  auto pref_int_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pref_int);
  int kbotbg_tmp = -1;
  for (int i = 0; i < static_cast<int>(pref_int_h.size()); ++i) {
    if (pref_int_h[i] < GWC::kbotbg_pref_max) {
      kbotbg_tmp = i;
    }
  }
  EKAT_REQUIRE_MSG(kbotbg_tmp >= 1,
    "Error! No interface pressure level found below kbotbg_pref_max. Check pref_int.\n");
  kbotbg_tmp -= 1; // move one level up

  s_common_init.do_molec_diff = do_molec_diff_in;
  s_common_init.nbot_molec = nbot_molec_in;
  s_common_init.ktop = ktop_in;
  s_common_init.kbotbg = kbotbg_tmp;
  s_common_init.kwv = kwv_in;
  s_common_init.oroko2 = GWC::half * kwv_in;
  s_common_init.effkwv = s_common_init.kwv * s_common_init.fcrit2;
  s_common_init.tndmax = s_common_init.orographic_only ? 500/GWC::sec_per_day : 400/GWC::sec_per_day;
}

} // namespace gw
} // namespace scream

#endif
