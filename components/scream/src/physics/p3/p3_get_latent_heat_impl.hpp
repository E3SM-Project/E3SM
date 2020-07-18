#ifndef P3_GET_LATENT_HEAT_IMPL_HPP
#define P3_GET_LATENT_HEAT_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
void Functions<S,D>
::get_latent_heat(const Int& ni, const Int& nk, view_2d<Spack>& v, view_2d<Spack>& s, view_2d<Spack>& f)
{
  using ExeSpace = typename KT::ExeSpace;

  constexpr Scalar lapvap  = C::LatVap;
  constexpr Scalar latice  = C::LatIce;

  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk);
  Kokkos::parallel_for("get_latent_heat", policy, KOKKOS_LAMBDA(const MemberType& team) {
    int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nk), [&] (const int& k) {
      v(i,k) = lapvap;
      s(i,k) = lapvap + latice;
      f(i,k) = latice;
    });
  });

  // Original P3 definition of latent heats:
  //       do i = its,ite
  //          do k = kts,kte
  //          xxlv(i,k)    = 3.1484e6-2370.*t(i,k)
  //          xxls(i,k)    = xxlv(i,k)+0.3337e6
  //          xlf(i,k)     = xxls(i,k)-xxlv(i,k)
  //          end do
  //       end do
}

} // namespace p3
} // namespace scream

#endif
