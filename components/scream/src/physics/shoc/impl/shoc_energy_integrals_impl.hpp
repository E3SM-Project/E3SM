#ifndef SHOC_ENERGY_INTEGRALS_IMPL_HPP
#define SHOC_ENERGY_INTEGRALS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_energy_integrals(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& host_dse,
  const uview_1d<const Spack>& pdel,
  const uview_1d<const Spack>& rtm,
  const uview_1d<const Spack>& rcm,
  const uview_1d<const Spack>& u_wind,
  const uview_1d<const Spack>& v_wind,
  Scalar&                      se_int,
  Scalar&                      ke_int,
  Scalar&                      wv_int,
  Scalar&                      wl_int)
{
  using ExeSpaceUtils = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  const auto ggr = C::gravit;

  // The team_barriers protect what we think is unexpected behavior in
  // Kokkos::parallel_reduce. We expect not to need these based on the semantics
  // of parallel_reduce. However, we speculate that in the Cuda implementation
  // of the team parallel_reduce, it's possible something in the bookkeeping at
  // the end of one parallel_reduce may sometimes be affecting the start of the
  // next, leading to nondeterminism. Then again, there might be something
  // subtly wrong in our code, and we're just not seeing it. In any case, these
  // team_barriers fix nondeterminism that was repeatedly and with high
  // confidence traced to these energy integral computations. We checked that
  // the view_reduction wrapper is not the cause by simplifying these to be bare
  // Kokkos::parallel_reduce calls acting on doubles and saw the same results.

  // Compute se_int
  se_int = ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {
    return host_dse(k)*pdel(k)/ggr;
  });
  team.team_barrier();

  // Compute ke_int
  ke_int = ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {
    return sp(0.5)*(ekat::square(u_wind(k))+ekat::square(v_wind(k)))*pdel(k)/ggr;
  });
  team.team_barrier();

  // Compute wv_int
  wv_int = ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {
    return (rtm(k)-rcm(k))*pdel(k)/ggr;
  });
  team.team_barrier();

  // Compute wl_int
  wl_int = ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {
    return rcm(k)*pdel(k)/ggr;
  });
  team.team_barrier();
}

} // namespace shoc
} // namespace scream

#endif
