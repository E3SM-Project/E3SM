#ifndef P3_FUNCTIONS_RAIN_SELF_COLLECTION_IMPL_HPP
#define P3_FUNCTIONS_RAIN_SELF_COLLECTION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_self_collection(const Spack& rho, const Spack& qr_incld, const Spack& nr_incld, Spack& nrslf)
{

}


} // namespace p3
} // namespace scream

#endif
