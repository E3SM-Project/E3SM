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
  const P3Runtime& runtime_options,
  const Smask& context)
{

  // Khroutdinov and Kogan (2000)
  const auto qc_not_small = qc_incld >= 1e-8 && context;

  const Scalar autoconversion_prefactor =
      runtime_options.autoconversion_prefactor;
  const Scalar autoconversion_qc_exponent =
      runtime_options.autoconversion_qc_exponent;
  const Scalar autoconversion_nc_exponent =
      runtime_options.autoconversion_nc_exponent;
  const Scalar autoconversion_radius = runtime_options.autoconversion_radius;

  Scalar CONS3;
  // TODO: always default to second branch after BFB stuff is addressed
  if(autoconversion_radius == sp(25.0e-6)) {
    CONS3 = C::CONS3;
  } else {
    CONS3 = sp(1.0) / (C::CONS2 * pow(autoconversion_radius, sp(3.0)));
  }

  if(qc_not_small.any()) {
    Spack sgs_var_coef;
    // sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, sp(2.47) );
    sgs_var_coef = 1;

    qc2qr_autoconv_tend.set(
        qc_not_small,
        sgs_var_coef * autoconversion_prefactor *
            pow(qc_incld, autoconversion_qc_exponent) *
            pow(nc_incld * sp(1.e-6) * rho, -autoconversion_nc_exponent));
    // note: ncautr is change in Nr; nc2nr_autoconv_tend is change in Nc
    ncautr.set(qc_not_small, qc2qr_autoconv_tend * CONS3);
    nc2nr_autoconv_tend.set(qc_not_small,
                            qc2qr_autoconv_tend * nc_incld / qc_incld);
  }

  nc2nr_autoconv_tend.set(qc2qr_autoconv_tend == 0 && context, 0);
  qc2qr_autoconv_tend.set(nc2nr_autoconv_tend == 0 && context, 0);
}

} // namespace p3
} // namespace scream

#endif
