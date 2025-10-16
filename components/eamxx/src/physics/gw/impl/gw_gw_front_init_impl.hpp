#ifndef GW_GW_FRONT_INIT_IMPL_HPP
#define GW_GW_FRONT_INIT_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_front_init. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::gw_front_init(
  // Inputs
  const Real& taubgnd,
  const Real& frontgfc_in,
  const Int& kfront_in)
{
  using exe_space_t = typename KT::ExeSpace;

  // The following are used to set the module data.

  // Parameters to calculate fav (average value of gaussian over bin).

  // Integration interval to get bin average.
  static constexpr Real dca  = 0.1;
  // Width of gaussian in phase speed.
  static constexpr Real c0   = 30;
  static constexpr Real half = 0.5;

  s_front_init.frontgfc = frontgfc_in;
  s_front_init.kfront = kfront_in;

  //! Allocate and calculate fav.
  const Int num_pgwv = s_common_init.pgwv*2 + 1;
  s_front_init.fav = view_1d<Real>("front.fav", num_pgwv);

  auto cref = s_common_init.cref;
  auto fav  = s_front_init.fav;
  auto dc   = s_common_init.dc;
  auto pgwv = s_common_init.pgwv;

  Kokkos::parallel_for(Kokkos::RangePolicy<exe_space_t>(0, num_pgwv), KOKKOS_LAMBDA(const Int l) {
    if (num_pgwv > 1) {
      //! Lower bound of bin.
      const Real cmn = cref(l) - half*dc;
      const Real cmx = cref(l) + half*dc;
      const Real cmnc0 = cmn/c0;
      const Real cmxc0 = cmx/c0;
      // Loop over integration intervals in bin.
      fav(l) = half * dca * (std::exp(-(cmnc0 * cmnc0)) + std::exp(-(cmxc0 * cmxc0)));
      const Int loop = static_cast<Int>(std::round(dc/dca));
      for (Int n = 1; n < loop; ++n) {
        const Real temp = (cmn+n*dca)/c0;
        fav(l) = fav(l) + dca * std::exp(-(temp*temp));
      }
      // Multiply by source strength.
      fav(l) = taubgnd * (fav(l)/dc);
    }

    // Prohibit wavenumber 0.
    fav(pgwv) = 0;
  });
}

} // namespace gw
} // namespace scream

#endif
