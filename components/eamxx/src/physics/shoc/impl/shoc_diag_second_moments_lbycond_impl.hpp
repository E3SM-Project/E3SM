#ifndef SHOC_DIAG_SECOND_MOMENTS_LBYCOND_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_LBYCOND_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_diag_second_moments_lbycond(
    const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc,
    const Scalar& ustar2, const Scalar& wstar,
    Scalar& wthl_sec, Scalar& wqw_sec, Scalar& uw_sec, Scalar& vw_sec,
    Scalar& wtke_sec, Scalar& thl_sec, Scalar& qw_sec, Scalar& qwthl_sec)
{
  // Purpose of this subroutine is to diagnose the lower
  // boundary condition for the second order moments needed
  // for the SHOC parameterization.
  // The thermodymnamic, tracer, and momentum fluxes are set
  // to the surface fluxes for the host model, while the
  // thermodynamic variances and covariances are computed
  // according to that of Andre et al. 1978.

  const Scalar ufmin = 0.01;

  auto uf = std::sqrt(ustar2+sp(0.3)*wstar*wstar);
  uf = ekat::impl::max<Scalar>(ufmin,uf);

  // Diagnose thermodynamics variances and covariances
  thl_sec   = sp(0.4)*sp(1.8)*((wthl_sfc/uf)*(wthl_sfc/uf));
  qw_sec    = sp(0.4)*sp(1.8)*((wqw_sfc/uf)*(wqw_sfc/uf)); 
  qwthl_sec = sp(0.2)*sp(1.8)*(wthl_sfc/uf)*(wqw_sfc/uf);

  // Vertical fluxes of heat and moisture, simply
  // use the surface fluxes given by host model
  wthl_sec = wthl_sfc;
  wqw_sec  = wqw_sfc;
  uw_sec   = uw_sfc;
  vw_sec   = vw_sfc;

  auto max_val = ekat::impl::max<Scalar>(sqrt(ustar2),ufmin);
  wtke_sec = max_val*max_val*max_val;
}
} // namespace shoc
} // namespace scream
#endif
