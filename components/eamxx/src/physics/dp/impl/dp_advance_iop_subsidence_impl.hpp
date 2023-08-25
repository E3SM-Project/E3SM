#ifndef DP_ADVANCE_IOP_SUBSIDENCE_IMPL_HPP
#define DP_ADVANCE_IOP_SUBSIDENCE_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace dp {

/*
 * Implementation of dp advance_iop_subsidence. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::do_advance_iop_subsidence_update(
    const Int& k,
    const Int& plev,
    const Spack& fac,
    const Spack& swfldint,
    const Spack& swfldint_p1,
    const uview_1d<const Spack>& in,
    const uview_1d<const Scalar>& in_s,
    const uview_1d<Spack>& update)
{
  Spack sin, sin_p1, sin_m1;

  auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
  auto range_pack2_m1_safe = range_pack1;
  auto range_pack2_p1_safe = range_pack1;
  range_pack2_m1_safe.set(range_pack1 < 1, 1); // don't want the shift to go below zero. we mask out that result anyway
  range_pack2_p1_safe.set(range_pack1 > plev-2, plev-2); // don't want the shift to go beyond pack.
  ekat::index_and_shift<-1>(in_s, range_pack2_m1_safe, sin, sin_m1);
  ekat::index_and_shift< 1>(in_s, range_pack2_p1_safe, sin, sin_p1);

  update(k) = in(k) - fac*(swfldint_p1*(sin_p1 - sin) + swfldint*(sin - sin_m1));
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::advance_iop_subsidence(
  const Int& plev,
  const Int& pcnst,
  const Scalar& scm_dt,
  const Scalar& ps_in,
  const uview_1d<const Spack>& u_in,
  const uview_1d<const Spack>& v_in,
  const uview_1d<const Spack>& t_in,
  const uview_2d<const Spack>& q_in,
  const uview_1d<const Spack>& hyai,
  const uview_1d<const Spack>& hyam,
  const uview_1d<const Spack>& hybi,
  const uview_1d<const Spack>& hybm,
  const uview_1d<const Spack>& wfld,
  const MemberType& team,
  const Workspace& workspace,
  const uview_1d<Spack>& u_update,
  const uview_1d<Spack>& v_update,
  const uview_1d<Spack>& t_update,
  const uview_2d<Spack>& q_update)
{
    // Local variables
  uview_1d<Spack>
    pmidm1, // pressure at model levels
    pintm1, // pressure at model interfaces (dim=plev+1)
    pdelm1, // pdel(k)   = pint  (k+1)-pint  (k)
    wfldint;// (dim=plev+1)
  workspace.template take_many_contiguous_unsafe<4>(
    {"pmidm1", "pintm1", "pdelm1", "wfldint"},
    {&pmidm1, &pintm1, &pdelm1, &wfldint});

  const Int plev_pack = ekat::npack<Spack>(plev);

  // Get vertical level profiles
  plevs0(plev, ps_in, hyai, hyam, hybi, hybm, team, pintm1, pmidm1, pdelm1);

  // Scalarize a bunch of views that need shift operations
  auto pmidm1_s  = scalarize(pmidm1);
  auto pintm1_s  = scalarize(pintm1);
  auto wfld_s    = scalarize(wfld);
  auto wfldint_s = scalarize(wfldint);
  auto u_in_s    = scalarize(u_in);
  auto v_in_s    = scalarize(v_in);
  auto t_in_s    = scalarize(t_in);

  wfldint_s(0) = 0;
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, plev_pack), [&] (Int k) {
      Spack spmidm1, spmidm1_m1, spintm1, spintm1_m1, swfld, swfld_m1;
      auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
      auto range_pack2 = range_pack1;
      range_pack2.set(range_pack1 < 1, 1); // don't want the shift to go below zero. we mask out that result anyway
      ekat::index_and_shift<-1>(pmidm1_s, range_pack2, spmidm1, spmidm1_m1);
      ekat::index_and_shift<-1>(pintm1_s, range_pack2, spintm1, spintm1_m1);
      ekat::index_and_shift<-1>(wfld_s,   range_pack2, swfld,   swfld_m1);
      Spack weight = (spintm1 - spintm1_m1) / (spmidm1 - spmidm1_m1);
      wfldint(k) = (1 - weight)*swfld_m1 + weight*swfld;
  });
  wfldint_s(plev) = 0;

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, plev_pack), [&] (Int k) {
      Spack swfldint, swfldint_p1;
      auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
      auto range_pack2 = range_pack1;
      range_pack2.set(range_pack1 > plev-1, plev-1);

      ekat::index_and_shift<1>(wfldint_s, range_pack2, swfldint, swfldint_p1);

      Spack fac = scm_dt/(2 * pdelm1(k));

      do_advance_iop_subsidence_update(k, plev, fac, swfldint, swfldint_p1, u_in, u_in_s, u_update);
      do_advance_iop_subsidence_update(k, plev, fac, swfldint, swfldint_p1, v_in, v_in_s, v_update);
      do_advance_iop_subsidence_update(k, plev, fac, swfldint, swfldint_p1, t_in, t_in_s, t_update);

      for (Int m = 0; m < pcnst; ++m) {
        // Grab m-th subview of q stuff
        auto q_update_sub = ekat::subview(q_update, m);
        auto q_in_sub     = ekat::subview(q_in, m);
        auto q_in_sub_s   = scalarize(q_in_sub);
        do_advance_iop_subsidence_update(k, plev, fac, swfldint, swfldint_p1, q_in_sub, q_in_sub_s, q_update_sub);
      }
  });

  // Top and bottom levels next
  Kokkos::Array<Int, 2> bot_top = {0, plev-1};
  for (Int i = 0; i < 2; ++i) {
    const auto k = bot_top[i];
    const auto pack_idx = ekat::npack<Spack>(k+1) - 1;
    const auto s_idx    = k % Spack::n;
    const auto idx1 = k == 0 ? k+1 : k;

    Scalar fac = scm_dt/(2 * pdelm1(pack_idx)[s_idx]);
    u_update(pack_idx)[s_idx] = u_in_s(k) - fac*(wfldint_s(idx1)*(u_in_s(idx1) - u_in_s(idx1-1)));
    v_update(pack_idx)[s_idx] = v_in_s(k) - fac*(wfldint_s(idx1)*(v_in_s(idx1) - v_in_s(idx1-1)));
    t_update(pack_idx)[s_idx] = t_in_s(k) - fac*(wfldint_s(idx1)*(t_in_s(idx1) - t_in_s(idx1-1)));

    for (Int m = 0; m < pcnst; ++m) {
      auto q_update_sub = ekat::subview(q_update, m);
      auto q_in_sub     = ekat::subview(q_in, m);
      auto q_in_sub_s   = scalarize(q_in_sub);

      q_update_sub(pack_idx)[s_idx] = q_in_sub_s(k) - fac*(wfldint_s(idx1)*(q_in_sub_s(idx1) - q_in_sub_s(idx1-1)));
    }
  }

  // thermal expansion term due to LS vertical advection
  constexpr Scalar rair = C::Rair;
  constexpr Scalar cpair = C::Cpair;
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, plev_pack), [&] (Int k) {
      t_update(k) = t_update(k) + scm_dt*wfld(k)*t_in(k)*rair/(cpair*pmidm1(k));
  });

  workspace.template release_many_contiguous<4>(
    {&pmidm1, &pintm1, &pdelm1, &wfldint});
}

} // namespace dp
} // namespace scream

#endif
