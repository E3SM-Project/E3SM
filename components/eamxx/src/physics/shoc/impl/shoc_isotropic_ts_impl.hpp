#ifndef SHOC_ISOTROPIC_TS_IMPL_HPP
#define SHOC_ISOTROPIC_TS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc isotropic_ts. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::isotropic_ts(
  const MemberType&            team,
  const Int&                   nlev,
  const Scalar&                lambda_low_in, 
  const Scalar&                lambda_high_in,
  const Scalar&                lambda_slope_in,
  const Scalar&                lambda_thresh_in,
  const Scalar&                brunt_int,
  const uview_1d<const Pack>& tke,
  const uview_1d<const Pack>& a_diss,
  const uview_1d<const Pack>& brunt,
  const uview_1d<Pack>&       isotropy)
{

  //constants from share/physics
  static constexpr  Scalar ggr = C::gravit.value;

  //Declare constants
         const     Scalar lambda_low    = lambda_low_in; 
         const     Scalar lambda_high   = lambda_high_in;
         const     Scalar lambda_slope  = lambda_slope_in;
         const     Scalar lambda_thresh = lambda_thresh_in;
  static constexpr Scalar maxiso        = 20000; // Return to isotropic timescale [s]

  const Int nlev_pack = ekat::npack<Pack>(nlev);

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

      // define the time scale
      const Pack tscale = 2*tke(k)/a_diss(k);

      // define a damping term "lambda" based on column stability
      Pack lambda(lambda_low + ((brunt_int/ggr)-lambda_thresh)*lambda_slope);
      lambda = ekat::max(Pack(lambda_low),ekat::min(Pack(lambda_high),lambda));

      const Pack buoy_sgs_save = brunt(k);
      lambda.set(buoy_sgs_save <=0, 0); //set lambda to zero where buoy_sgs_save <=0

      // Compute the return to isotropic timescale
      isotropy(k)= ekat::min(Pack(maxiso),tscale/(1+lambda*buoy_sgs_save*ekat::square(tscale)));

    });
}

} // namespace shoc
} // namespace scream

#endif
