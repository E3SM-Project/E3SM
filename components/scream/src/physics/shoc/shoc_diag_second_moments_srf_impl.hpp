
#ifndef SHOC_DIAG_SECOND_MOMENTS_SRF_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_SRF_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_diag_second_moments_srf(
    const Int& shcol, 
    const Spack& wthl_sfc, const Spack& uw_sfc, const Spack& vw_sfc,
    Spack& ustar2, Spack& wstar)
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

 ustar2 = pack::sqrt(uw_sfc*uw_sfc+vw_sfc*vw_sfc);

 const auto is_wthl_ge_zero = wthl_sfc >= zero;

 wstar.set(is_wthl_ge_zero,
           pack::pow(one/basetemp*ggr*wthl_sfc*one, third));

 wstar.set(!is_wthl_ge_zero, zero);

}

} // namespace shoc
} // namespace scream

#endif
