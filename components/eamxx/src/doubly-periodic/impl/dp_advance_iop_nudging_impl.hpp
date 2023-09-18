#ifndef DP_ADVANCE_IOP_NUDGING_IMPL_HPP
#define DP_ADVANCE_IOP_NUDGING_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp advance_iop_nudging. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::advance_iop_nudging(
    // Input arguments
    const Int& plev,
    const Scalar& scm_dt,
    const Scalar& ps_in,
    const uview_1d<const Spack>& t_in,
    const uview_1d<const Spack>& q_in,
    const uview_1d<const Spack>& tobs,
    const uview_1d<const Spack>& qobs,
    const uview_1d<const Spack>& hyai,
    const uview_1d<const Spack>& hyam,
    const uview_1d<const Spack>& hybi,
    const uview_1d<const Spack>& hybm,
    // Kokkos stuff
    const MemberType& team,
    const Workspace& workspace,
    // Output arguments
    const uview_1d<Spack>& t_update,
    const uview_1d<Spack>& q_update,
    const uview_1d<Spack>& relaxt,
    const uview_1d<Spack>& relaxq)
{
  // Local variables
  uview_1d<Spack>
    pmidm1, // pressure at model levels
    pintm1, // pressure at model interfaces (dim=plev+1)
    pdelm1, // pdel(k)   = pint  (k+1)-pint  (k)
    rtau;
  workspace.template take_many_contiguous_unsafe<4>(
    {"pmidm1", "pintm1", "pdelm1", "rtau"},
    {&pmidm1, &pintm1, &pdelm1, &rtau});

  // Get vertical level profiles
  plevs0(plev, ps_in, hyai, hyam, hybi, hybm, team, pintm1, pmidm1, pdelm1);

  const Int plev_pack = ekat::npack<Spack>(plev);
  constexpr Scalar iop_nudge_tq_low = DPC::iop_nudge_tq_low;
  constexpr Scalar iop_nudge_tq_high = DPC::iop_nudge_tq_high;
  constexpr Scalar iop_nudge_tscale = DPC::iop_nudge_tscale;


  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, plev_pack), [&] (Int k) {
    relaxt(k) = 0;
    relaxq(k) = 0;

    const auto condition = pmidm1(k) <= iop_nudge_tq_low*100
                           &&
                           pmidm1(k) >= iop_nudge_tq_high*100;
    rtau(k).set(condition, ekat::impl::max(scm_dt, iop_nudge_tscale));
    relaxt(k).set(condition, (t_in(k) - tobs(k))/rtau(k));
    relaxq(k).set(condition, (q_in(k) - qobs(k))/rtau(k));

    t_update(k).set(condition, t_in(k) + relaxt(k)*scm_dt);
    q_update(k).set(condition, q_in(k) + relaxq(k)*scm_dt);
  });
}

} // namespace dp
} // namespace scream

#endif
