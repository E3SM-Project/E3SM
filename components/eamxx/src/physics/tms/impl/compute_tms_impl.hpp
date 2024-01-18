#ifndef COMPUTE_TMS_IMPL_HPP
#define COMPUTE_TMS_IMPL_HPP

#include "tms_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace tms {

template<typename S, typename D>
void Functions<S,D>::compute_tms(
  const int&                   ncols,
  const int&                   nlevs,
  const view_3d<const Scalar>& horiz_wind,
  const view_2d<const Scalar>& t_mid,
  const view_2d<const Scalar>& p_mid,
  const view_2d<const Scalar>& exner,
  const view_2d<const Scalar>& z_mid,
  const view_1d<const Scalar>& sgh,
  const view_1d<const Scalar>& landfrac,
  const view_1d<Scalar>&       ksrf,
  const view_2d<Scalar>&       tau_tms)
{
  using C  = physics::Constants<Scalar>;

  // Define some constants used
  const Scalar horomin = 1;          // Minimum value of subgrid orographic height for mountain stress [ m ]
  const Scalar z0max   = 100;        // Maximum value of z_0 for orography [ m ]
  const Scalar dv2min  = 0.01;       // Minimum shear squared [ m2/s2 ]
  const Scalar orocnst = C::orocnst; // Converts from standard deviation to height [ no unit ]
  const Scalar z0fac   = C::z0fac;   // Factor determining z_0 from orographic standard deviation [ no unit ]
  const Scalar karman  = C::Karman;  // von Karman constant
  const Scalar gravit  = C::gravit;  // Acceleration due to gravity
  const Scalar rair    = C::Rair;    // Gas constant for dry air

  // Loop over columns
  const typename KT::RangePolicy policy (0,ncols);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int& i) {
    // Subview on column, scalarize since we only care about last 2 levels (never loop over levels)
    const auto u_wind_i = ekat::subview(horiz_wind, i, 0);
    const auto v_wind_i = ekat::subview(horiz_wind, i, 1);
    const auto t_mid_i  = ekat::subview(t_mid, i);
    const auto p_mid_i  = ekat::subview(p_mid, i);
    const auto exner_i  = ekat::subview(exner, i);
    const auto z_mid_i  = ekat::subview(z_mid, i);

    // Determine subgrid orgraphic height (mean to peak)
    const auto horo = orocnst*sgh(i);

    if (horo < horomin) {
      // No mountain stress if horo is too small
      ksrf(i)       = 0;
      tau_tms(i, 0) = 0;
      tau_tms(i, 1) = 0;
    } else {
      // Determine z0m for orography
      const auto z0oro = ekat::impl::min(z0fac*horo, z0max);

      // Calculate neutral drag coefficient
      const auto tmp = karman/std::log((z_mid_i(nlevs-1) + z0oro )/z0oro);
      auto cd = tmp*tmp;

      // Calculate the Richardson number over the lowest 2 layers
      const auto kb = nlevs-1;
      const auto kt = nlevs-2;
      const auto tmp_u = u_wind_i(kt) - u_wind_i(kb);
      const auto tmp_v = v_wind_i(kt) - v_wind_i(kb);
      const auto dv2 = ekat::impl::max(tmp_u*tmp_u + tmp_v*tmp_v, dv2min);

      const auto ri  = 2*gravit*(t_mid_i(kt)*exner_i(kt) - t_mid_i(kb)*exner_i(kb))*(z_mid_i(kt) - z_mid_i(kb))/
                      ((t_mid_i(kt)*exner_i(kt) + t_mid_i(kb)*exner_i(kb))*dv2);

      // Calculate the instability function and modify the neutral drag cofficient.
      // We should probably follow more elegant approach like Louis et al (1982) or
      // Bretherton and Park (2009) but for now we use very crude approach : just 1
      // for ri < 0, 0 for ri > 1, and linear ramping.
      const auto stabfri = ekat::impl::max(0.0, ekat::impl::min(1.0, 1.0-ri));
      cd *= stabfri;

      // Compute density, velocity magnitude and stress using bottom level properties
      const auto rho = p_mid_i(nlevs-1)/(rair*t_mid_i(nlevs-1));
      const auto vmag = std::sqrt(u_wind_i(nlevs-1)*u_wind_i(nlevs-1) +
                                  v_wind_i(nlevs-1)*v_wind_i(nlevs-1));
      ksrf(i)       = rho*cd*vmag*landfrac(i);
      tau_tms(i, 0) = -ksrf(i)*u_wind_i(nlevs-1);
      tau_tms(i, 1) = -ksrf(i)*v_wind_i(nlevs-1);
    }
  });
}

} // namespace scream
} // namespace tms

#endif // COMPUTE_TMS_IMPL_HPP
