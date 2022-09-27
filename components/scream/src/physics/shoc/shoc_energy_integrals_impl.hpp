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

  Scalar se_tmp = 0, ke_tmp = 0, wv_tmp = 0, wl_tmp = 0;
  const auto shost_dse = ekat::scalarize(host_dse);
  const auto spdel     = ekat::scalarize(pdel);
  const auto srtm      = ekat::scalarize(rtm);
  const auto srcm      = ekat::scalarize(rcm);
  const auto su_wind   = ekat::scalarize(u_wind);

  const auto sv_wind   = ekat::scalarize(v_wind);
  // Compute se_int
  ExeSpaceUtils::parallel_reduce(team,0,nlev,
                                 [&] (const int k, Scalar& local_sum) {
                                   local_sum += shost_dse(k)*spdel(k)/ggr;
                                 }, se_tmp);
  team.team_barrier();

  // Compute ke_int
  ExeSpaceUtils::parallel_reduce(team,0,nlev,
                                 [&] (const int k, Scalar& local_sum) {
                                   local_sum += sp(0.5) * ((su_wind(k)*su_wind(k)) + (sv_wind(k)*sv_wind(k))) * spdel(k)/ggr;
                                 }, ke_tmp);
  team.team_barrier();

  // Compute wv_int
  ExeSpaceUtils::parallel_reduce(team,0,nlev,
                                 [&] (const int k, Scalar& local_sum) {
                                   local_sum += (srtm(k)-srcm(k))*spdel(k)/ggr;
                                 }, wv_tmp);
  team.team_barrier();

  // Compute wl_int
  ExeSpaceUtils::parallel_reduce(team,0,nlev,
                                 [&] (const int k, Scalar& local_sum) {
                                   local_sum += srcm(k)*spdel(k)/ggr;
                                 }, wl_tmp);
  team.team_barrier();

  Kokkos::single(Kokkos::PerTeam(team), [&] () {
      se_int = se_tmp;
      ke_int = ke_tmp;
      wv_int = wv_tmp;
      wl_int = wl_tmp;
  });
}

} // namespace shoc
} // namespace scream

#endif
