#ifndef SHOC_DIAG_SECOND_MOMENTS_UBYCOND_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_UBYCOND_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_diag_second_moments_ubycond(
    Scalar& thl_sec, Scalar& qw_sec, Scalar& wthl_sec, Scalar& wqw_sec,
    Scalar& qwthl_sec, Scalar& uw_sec, Scalar& vw_sec, Scalar& wtke_sec)
{
  // Purpose of this subroutine is to diagnose the upper
  //  boundary condition for the second order moments
  //  needed for the SHOC parameterization.  Currently
  //  set all to zero.

 const auto zero = C::ZERO;

 wthl_sec  = zero;
 wqw_sec   = zero;
 uw_sec    = zero;
 vw_sec    = zero;
 wtke_sec  = zero;
 thl_sec   = zero;
 qw_sec    = zero;
 qwthl_sec = zero;

}

} // namespace shoc
} // namespace scream
#endif
