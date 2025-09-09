#ifndef GW_GW_EDIFF_IMPL_HPP
#define GW_GW_EDIFF_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "util/eamxx_utils.hpp"
#include <ekat_math_utils.hpp>
#include <ekat_subview_utils.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_ediff. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_ediff(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const Int& pver,
  const Int& pgwv,
  const Int& kbot,
  const Int& ktop,
  const Int& tend_level,
  const Real& dt,
  const uview_2d<const Real>& gwut,
  const uview_1d<const Real>& ubm,
  const uview_1d<const Real>& nm,
  const uview_1d<const Real>& rho,
  const uview_1d<const Real>& pmid,
  const uview_1d<const Real>& rdpm,
  const uview_1d<const Real>& c,
  // Outputs
  const uview_1d<Real>& egwdffi,
  const uview_1d<Real>& decomp_ca,
  const uview_1d<Real>& decomp_cc,
  const uview_1d<Real>& decomp_dnom,
  const uview_1d<Real>& decomp_ze)
{
  //-----------------------------Local Workspace------------------------------

  // 0._r8 (array to pass to vd_lu_decomp).
  Real zero = 0;
  // Effective gw diffusivity at midpoints.
  uview_1d<Real> egwdffm, tmpi2;
  workspace.template take_many_contiguous_unsafe<2>(
    {"egwdffm", "tmpi2"},
    {&egwdffm, &tmpi2});

  // Inverse Prandtl number.
  static constexpr Real prndl=0.25;
  // Density scale height.
  static constexpr Real dscale=7000;
  static constexpr Real half=0.5;

  //--------------------------------------------------------------------------

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver+1), [&] (const int k) {
      egwdffi(k) = 0;
      egwdffm(k) = 0;
      tmpi2(k)   = 0;
  });

  // Calculate effective diffusivity at midpoints.
  const int num_pgwv = 2*pgwv + 1;
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, num_pgwv * (ktop+1), num_pgwv * (kbot+1)), [&] (const int k_pgwv) {
      const int k = k_pgwv / num_pgwv;
      const int l = k_pgwv % num_pgwv;

      const Real add = prndl * half * gwut(k,l) * (c(l)-ubm(k)) / bfb_square(nm(k));
      Kokkos::atomic_add(&egwdffm(k), add);
  });

  // Interpolate effective diffusivity to interfaces.
  // Assume zero at top and bottom interfaces.
  midpoint_interp(team,
                  ekat::subview(egwdffm, Kokkos::pair<int, int>{ktop+1, kbot+1}),
                  ekat::subview(egwdffi, Kokkos::pair<int, int>{ktop+2, kbot+1}));

  // Limit diffusivity to some reasonable value.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver+1), [&] (const int k) {
      egwdffi(k) = ekat::impl::min((Real)150., egwdffi(k));
  });

  // Do not calculate diffusivities below level where tendencies are
  // actually allowed.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, tend_level+1, kbot+1), [&] (const int k) {
      egwdffi(k) = 0;
  });

  // Calculate dt * (gravit*rho)^2/dp at interior interfaces.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, ktop+2, kbot+2), [&] (const int k) {
      tmpi2(k) = dt * bfb_square(C::gravit * rho(k)) / (pmid(k) - pmid(k-1));
  });

  // Decompose the diffusion matrix.
  // Note that [ktop,kbot] are model interfaces (beginning at zero), whereas
  // in vd_lu_decomp they are expected as midpoints.
  vd_lu_decomp(team, pver,
               zero, egwdffi, tmpi2, rdpm, dt, zero, ktop+1, kbot+1, decomp_ca, decomp_cc, decomp_dnom, decomp_ze);

  workspace.template release_many_contiguous<2>(
    {&egwdffm, &tmpi2});
}

} // namespace gw
} // namespace scream

#endif
