#ifndef P3_SUBGRID_VARIANCE_SCALING_IMPL_HPP
#define P3_SUBGRID_VARIANCE_SCALING_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::subgrid_variance_scaling(const Spack& relvar, const Scalar& expon)
{
  /* We assume subgrid variations in qc follow a gamma distribution with inverse
     relative variance relvar = 1/(var(qc)/qc**2). In this case, if the tendency
     for a given process is of the form A*qc**expon for a local value of qc, then
     the average process rate over the PDF is
     gamma(relvar+expon)/[gamma(relvar)*relvar**expon]*A*average(qc)**expon. This
     function calculates the local process rate => cell-average process rate scaling
     factor gamma(relvar+expon)/[gamma(relvar)*relvar**expon]. See Morrison and
     Gettelman (2008; JCLI) eq 9 for details.
  */

  /* ***bounds checking not operational yet ***

  // Check that relvar is within allowable range
  //============================================
  const Scalar relvar_min = 0.1;
  const Scalar relvar_max = 10.0;

  const auto relvar_exceeds_bounds  = !(relvar > relvar_min && relvar < relvar_max);

  if (relvar_exceeds_bounds.any()) {
    EKAT_REQUIRE_MSG( condition, "relvar outside allowable bounds" );
  }

  // Check that expon >0.
  //============================================
  if (expon < 0.0){
    const auto msg = "expon<0. This might be ok, but isn't unit tested and can drive subgrid_variance_scaling negative. Be careful if you proceed.";
    EKAT_REQUIRE_MSG( condition, msg );
  }

  */

  // Compute result
  //============================================
  Spack result;
  Spack exponent=Spack(expon);
  result=tgamma( relvar+exponent ) / ( tgamma(relvar)*pow(relvar,exponent) );

  return result;

}

} // namespace p3
} // namespace scream

#endif // P3_SUBGRID_VARIANCE_SCALING_IMPL_HPP
