#ifndef P3_FUNCTIONS_CLOUD_SED_IMPL_HPP
#define P3_FUNCTIONS_CLOUD_SED_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 cloud sedimentation function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_sedimentation(
    const uview_1d<const Spack>& qc_incld,
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& lcldm,
    const uview_1d<const Spack>& acn,
    const uview_1d<const Spack>& inv_dzq,
    const MemberType& team,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& odt, const bool& log_predictNc,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& mu_c,
    const uview_1d<Spack>& lamc,
    const uview_1d<Spack>& qc_tend,
    const uview_1d<Spack>& nc_tend,
    Scalar& prt_liq)
{
  // TODO
}

} // namespace p3
} // namespace scream

#endif
