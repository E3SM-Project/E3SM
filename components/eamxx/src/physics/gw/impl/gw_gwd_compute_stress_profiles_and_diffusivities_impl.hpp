#ifndef GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP
#define GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_compute_stress_profiles_and_diffusivities. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_compute_stress_profiles_and_diffusivities(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const Int& ngwv,
const uview_1d<const Int>& src_level,
const uview_1d<const Spack>& ubi,
const uview_1d<const Spack>& c,
const uview_1d<const Spack>& rhoi,
const uview_1d<const Spack>& ni,
const uview_1d<const Spack>& kvtt,
const uview_1d<const Spack>& t,
const uview_1d<const Spack>& ti,
const uview_1d<const Spack>& piln,
// Inputs/Outputs
const uview_1d<Spack>& tau)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
