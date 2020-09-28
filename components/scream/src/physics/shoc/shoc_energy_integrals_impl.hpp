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

  // Use serialized version of view_reduction below if doing BFB testing
#ifdef EKAT_DEFAULT_BFB
  constexpr bool Serialize = true;
#else
  constexpr bool Serialize = false;
#endif

  // Compute se_int
  ExeSpaceUtils::template view_reduction<Serialize>(team,0,nlev,
                                                    [&] (const int k) -> Spack {
    return host_dse(k)*pdel(k)/ggr;
  }, se_int);

  // Compute ke_int
  ExeSpaceUtils::template view_reduction<Serialize>(team,0,nlev,
                                                    [&] (const int k) -> Spack {
    return sp(0.5)*(ekat::square(u_wind(k))+ekat::square(v_wind(k)))*pdel(k)/ggr;
  }, ke_int);

  // Compute wv_int
  ExeSpaceUtils::template view_reduction<Serialize>(team,0,nlev,
                                                    [&] (const int k) -> Spack {
    return (rtm(k)-rcm(k))*pdel(k)/ggr;
  }, wv_int);

  // Compute wl_int
  ExeSpaceUtils::template view_reduction<Serialize>(team,0,nlev,
                                                    [&] (const int k) -> Spack {
    return rcm(k)*pdel(k)/ggr;
  }, wl_int);
}

} // namespace shoc
} // namespace scream

#endif
