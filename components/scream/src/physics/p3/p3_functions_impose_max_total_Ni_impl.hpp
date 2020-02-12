#ifndef P3_FUNCTIONS_AUTOCONVERSION_IMPL_HPP
#define P3_FUNCTIONS_AUTOCONVERSION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::impose_max_total_Ni(Spack& nitot_local, Spack& max_total_Ni, Spack& inv_rho_local)
{
    //--------------------------------------------------------------------------------
    // Impose maximum total ice number concentration (total of all ice categories).
    // If the sum of all nitot(:) exceeds maximum allowable, each category to preserve
    // ratio of number between categories.
    //--------------------------------------------------------------------------------


}

} // namespace p3
} // namespace scream



#endif