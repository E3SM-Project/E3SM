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
  Kokkos::View<Real*, Kokkos::HostSpace> pref_int,
  const bool& do_molec_diff_in,
  const Int& nbot_molec_in,
  const Int& ktop_in,
  const Real& kwv_in)
{
  EKAT_REQUIRE_MSG(pref_int.size()==pver_in+1, "Error! pref_int size is incorrect.\n");

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

  // // calculate reference pressure at interfaces
  // const auto hyai = m_grid->get_geometry_data("hyai").get_view<const Real*>();
  // const auto hybi = m_grid->get_geometry_data("hybi").get_view<const Real*>();
  // Kokkos::View<Real*, Kokkos::HostSpace> pref_int("pref_int", hyai.size());
  // // Kokkos::parallel_for(Kokkos::RangePolicy<>(0, m_nlev+1), KOKKOS_LAMBDA (const int k) {
  // Kokkos::parallel_for("gw_common_init_calculate_pref_int", 
  //   Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, hyai.size()), 
  //   KOKKOS_LAMBDA (const int k) {
  //     pref_int(k) = PC::P0.value*hyai(k) + PC::P0.value*hybi(k);
  // });
  // Kokkos::fence();

  // Set phase speeds
  for (int l = -s_common_init.pgwv; l <= s_common_init.pgwv; ++l)
    s_common_init.cref[l + s_common_init.pgwv] = s_common_init.dc * l;

  // cref is a device view; initialize via a HostSpace mirror, then deep_copy.
  auto cref_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), s_common_init.cref);
  for (int l = -s_common_init.pgwv; l <= s_common_init.pgwv; ++l) {
    cref_h[l + s_common_init.pgwv] = s_common_init.dc * l;
  }
  Kokkos::deep_copy(s_common_init.cref, cref_h);

  s_common_init.orographic_only = !s_common_init.use_gw_convect && !s_common_init.use_gw_frontal;

  // calculate alpha (Newtonian cooling coefficients)
  s_common_init.alpha = view_1d<Real>("alpha", GWC::nalph);
  if (s_common_init.orographic_only) {
    Kokkos::deep_copy(s_common_init.alpha,1e-6);
  } else {
    // pre-calculated newtonian damping:
    //     * convert alpha0 from 1/day to 1/s
    //     * ensure alpha0 is not smaller than 1e-6
    //     * convert alpha_pressure_mb to pa
    // Create Host Views
    Kokkos::View<Real*, Kokkos::HostSpace> alpha0_per_sec("alpha0_per_sec", GWC::nalph);
    Kokkos::View<Real*, Kokkos::HostSpace> alpha_pressure_pa("alpha_pressure_pa", GWC::nalph);
    // Newtonian damping calculation in a parallel kernel (on Host)
    Kokkos::parallel_for("gw_common_init_calculate_newtonian_damping", 
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, GWC::nalph), 
      KOKKOS_LAMBDA (const int k) {
        alpha0_per_sec(k) = std::max(GWC::alpha0[k] / 86400.0, 1.0e-6);
        alpha_pressure_pa(k) = GWC::alpha_pressure_mb[k] * 1.0e2;
      }
    );
    // Interpolate directly into the Host View
    Kokkos::View<Real*, Kokkos::HostSpace> alpha_tmp("alpha_tmp", pver_in+1);
    simple_linear_interp(alpha_pressure_pa, alpha0_per_sec, pref_int, alpha_tmp);
    // copy to device struct
    Kokkos::deep_copy(s_common_init.alpha, alpha_tmp);
  }

  // set bottom index for background spectrum
  int kbotbg_tmp = -1;
  for (int i = 0; i < pref_int.size(); ++i) {
    if (pref_int[i] < GWC::kbotbg_pref_max) {
      kbotbg_tmp = i;
    }
  }
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
