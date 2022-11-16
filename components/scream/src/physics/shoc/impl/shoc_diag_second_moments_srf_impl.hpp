
#ifndef SHOC_DIAG_SECOND_MOMENTS_SRF_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_SRF_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_diag_second_moments_srf(
    const Scalar& wthl_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc,
    Scalar& ustar2, Scalar& wstar)
{
 // Purpose of this subroutine is to diagnose surface
 // properties needed for the the lower
 // boundary condition for the second order moments needed
 // for the SHOC parameterization.

 const auto one      = C::ONE;
 const auto zero     = C::ZERO;
 const auto third    = C::THIRD;
 const auto ggr      = C::gravit;
 const auto basetemp = C::basetemp;

 ustar2 = std::sqrt(uw_sfc*uw_sfc+vw_sfc*vw_sfc);

 const auto is_wthl_ge_zero = wthl_sfc >= zero;

 if (is_wthl_ge_zero) {
   wstar = std::pow(one/basetemp*ggr*wthl_sfc*one, third);
 }
 else {
   wstar = zero;
 }
}
} // namespace shoc
} // namespace scream

#endif
