#ifndef P3_FUNCTIONS_EVAPORATE_SUBLIMATE_PRECIP_IMPL.HPP
#define P3_FUNCTIONS_EVAPORATE_SUBLIMATE_PRECIP_IMPL.HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::evaporate_sublimate_precip(const Spack& qr_incld, const Spack& qc_incld, const Spack& nr_incld, const Spack& qitot_incld,
			     const Spack& lcldm, const Spack& rcldm, const Spack& qvs, const Spack& ab, const Spack& epsr,
			     const Spack& qv, Spack& qrevp, Spack& nrevp)
{
}


} // namespace p3
} // namespace scream

#endif

