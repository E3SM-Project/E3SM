#ifndef P3_ICE_SUPERSAT_CONSERVATION_IMPL_HPP
#define P3_ICE_SUPERSAT_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice_supersat_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::ice_supersat_conservation(Spack& qidep, Spack& qinuc, const Spack& cld_frac_i, const Spack& qv, const Spack& qv_sat_i, const Spack& latent_heat_sublim, const Spack& t_atm, const Real& dt)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace p3
} // namespace scream

#endif
