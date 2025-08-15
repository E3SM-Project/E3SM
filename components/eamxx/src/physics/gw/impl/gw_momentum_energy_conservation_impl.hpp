#ifndef GW_MOMENTUM_ENERGY_CONSERVATION_IMPL_HPP
#define GW_MOMENTUM_ENERGY_CONSERVATION_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw momentum_energy_conservation. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::momentum_energy_conservation(
// Inputs
const Int& pver,
const Int& ncol,
const uview_1d<const Int>& tend_level,
const Spack& dt,
const uview_1d<const Spack>& taucd,
const uview_1d<const Spack>& pint,
const uview_1d<const Spack>& pdel,
const uview_1d<const Spack>& u,
const uview_1d<const Spack>& v,
// Inputs/Outputs
const uview_1d<Spack>& dudt,
const uview_1d<Spack>& dvdt,
const uview_1d<Spack>& dsdt,
const uview_1d<Spack>& utgw,
const uview_1d<Spack>& vtgw,
const uview_1d<Spack>& ttgw)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
