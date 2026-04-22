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
      const Real cmn = cref(l) - GWC::half*dc;
      const Real cmx = cref(l) + GWC::half*dc;
      const Real cmnc0 = cmn/GWC::c0;
      const Real cmxc0 = cmx/GWC::c0;
      // Loop over integration intervals in bin.
      fav(l) = GWC::half * GWC::dca * (std::exp(-(cmnc0 * cmnc0)) + std::exp(-(cmxc0 * cmxc0)));
      const Int loop = static_cast<Int>(std::round(dc/GWC::dca));
      for (Int n = 1; n < loop; ++n) {
        const Real temp = (cmn+n*GWC::dca)/GWC::c0;
        fav(l) = fav(l) + GWC::dca * std::exp(-(temp*temp));
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
