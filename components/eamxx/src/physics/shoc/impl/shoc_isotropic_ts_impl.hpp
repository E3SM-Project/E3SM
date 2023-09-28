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
  const Scalar&                brunt_int,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& a_diss,
  const uview_1d<const Spack>& brunt,
  const uview_1d<Spack>&       isotropy)
{

  //constants from physics/share
  static constexpr  Scalar ggr = C::gravit;

  //Declare constants
  static constexpr Scalar lambda_low   = 0.001;
  static constexpr Scalar lambda_high  = 0.04;
  static constexpr Scalar lambda_slope = 2.65;
  static constexpr Scalar lambda_thresh= 0.02;
  static constexpr Scalar maxiso       = 20000; // Return to isotropic timescale [s]

  const Int nlev_pack = ekat::npack<Spack>(nlev);

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

      // define the time scale
      const Spack tscale = 2*tke(k)/a_diss(k);

      // define a damping term "lambda" based on column stability
      Spack lambda(lambda_low + ((brunt_int/ggr)-lambda_thresh)*lambda_slope);
      lambda = ekat::max(Spack(lambda_low),ekat::min(Spack(lambda_high),lambda));

      const Spack buoy_sgs_save = brunt(k);
      lambda.set(buoy_sgs_save <=0, 0); //set lambda to zero where buoy_sgs_save <=0

      // Compute the return to isotropic timescale
      isotropy(k)= ekat::min(Spack(maxiso),tscale/(1+lambda*buoy_sgs_save*ekat::square(tscale)));

    });
}

} // namespace shoc
} // namespace scream

#endif
