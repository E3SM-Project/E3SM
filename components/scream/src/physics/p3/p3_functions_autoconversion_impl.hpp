#ifndef P3_FUNCTIONS_AUTOCONVERSION_IMPL_HPP
#define P3_FUNCTIONS_AUTOCONVERSION_IMPL_HPP

#include "p3_functions.hpp"
#include "p3_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_water_autoconversion(const Spack& rho, const Spack& qc_incld, const Spack& nc_incld,
    Spack& qcaut, Spack& ncautc, Spack& ncautr)
{
    // Khroutdinov and Kogan (2000)
    const auto qc_not_small = qc_incld >=  1e-8;
    constexpr Scalar CONS3 = C::CONS3;
    if(qc_not_small.any()){
        qcaut.set(qc_not_small,
            sp(1350.0)*pow(qc_incld,sp(2.47))*pow(nc_incld*sp(1.e-6)*rho,sp(-1.79)));
        // note: ncautr is change in Nr; ncautc is change in Nc
        ncautr.set(qc_not_small, qcaut*CONS3);
        ncautc.set(qc_not_small, qcaut*nc_incld/qc_incld);
    }

    ncautc.set(qcaut == 0, 0.0);
    qcaut.set(ncautc == 0, 0.0);

}

} // namespace p3
} // namespace scream

#endif
