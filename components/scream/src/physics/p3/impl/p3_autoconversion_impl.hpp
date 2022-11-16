#ifndef P3_AUTOCONVERSION_IMPL_HPP
#define P3_AUTOCONVERSION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_subgrid_variance_scaling_impl.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_water_autoconversion(
  const Spack& rho, const Spack& qc_incld, const Spack& nc_incld,
  const Spack& inv_qc_relvar, Spack& qc2qr_autoconv_tend, Spack& nc2nr_autoconv_tend, Spack& ncautr,
  const Smask& context)
{
  // Khroutdinov and Kogan (2000)
  const auto qc_not_small = qc_incld >= 1e-8 && context;
  constexpr Scalar CONS3 = C::CONS3;
  if(qc_not_small.any()){
    Spack sgs_var_coef;
    // sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, sp(2.47) );
    sgs_var_coef = 1;

    qc2qr_autoconv_tend.set(qc_not_small,
              sgs_var_coef*1350*pow(qc_incld,sp(2.47))*pow(nc_incld*sp(1.e-6)*rho,sp(-1.79)));
    // note: ncautr is change in Nr; nc2nr_autoconv_tend is change in Nc
    ncautr.set(qc_not_small, qc2qr_autoconv_tend*CONS3);
    nc2nr_autoconv_tend.set(qc_not_small, qc2qr_autoconv_tend*nc_incld/qc_incld);
  }

  nc2nr_autoconv_tend.set(qc2qr_autoconv_tend == 0 && context, 0);
  qc2qr_autoconv_tend.set(nc2nr_autoconv_tend == 0 && context, 0);
}

} // namespace p3
} // namespace scream

#endif
