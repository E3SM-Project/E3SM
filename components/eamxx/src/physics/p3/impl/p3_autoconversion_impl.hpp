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
  const Spack& rho, const Spack& qc_incld, const Spack& qr_incld, const Spack& nc_incld,
  const Spack& inv_qc_relvar, Spack& qc2qr_autoconv_tend, Spack& nc2nr_autoconv_tend, Spack& ncautr, Spack& nu, Spack& inv_rho,
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

  const Scalar CONS3 = sp(1.0) / (C::CONS2 * pow(autoconversion_radius, sp(3.0)));
 
  constexpr Scalar kc = C::kc;

  const bool use_KK = runtime_options.use_KK;

  if(qc_not_small.any()) {
    Spack sgs_var_coef;
    // sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, sp(2.47) );
    sgs_var_coef = 1;
    if(use_KK){
      qc2qr_autoconv_tend.set(
        qc_not_small,
        sgs_var_coef * autoconversion_prefactor *
            pow(qc_incld, autoconversion_qc_exponent) *
            pow(nc_incld * sp(1.e-6) * rho, -autoconversion_nc_exponent));
    // note: ncautr is change in Nr; nc2nr_autoconv_tend is change in Nc
      ncautr.set(qc_not_small, qc2qr_autoconv_tend * CONS3);
      nc2nr_autoconv_tend.set(qc_not_small,
                            qc2qr_autoconv_tend * nc_incld / qc_incld);}
   else { //SB (2003)
      Spack dum;
      Spack dum1;
      dum   = 1.0 - (qc_incld/(qc_incld+qr_incld));
      dum1  = 600.*pow(dum,0.68)*pow(1.-pow(dum,0.68),3);
      qc2qr_autoconv_tend.set(qc_not_small,kc*1.9230769e-5*(nu+2.)*(nu+4.)/pow(nu+1.,2)*         
                      pow(rho*qc_incld*1.e-3,4)/pow(rho*nc_incld*1.e-6,2)*(1.+dum1/pow(1.-dum,2))*1000.*inv_rho);
    // note: ncautr is change in Nr; nc2nr_autoconv_tend is change in Nc
      ncautr.set(qc_not_small, qc2qr_autoconv_tend * 7.6923076e+9);
      nc2nr_autoconv_tend.set(qc_not_small,
                            qc2qr_autoconv_tend * 0.5);}

  }
  nc2nr_autoconv_tend.set(qc2qr_autoconv_tend == 0 && context, 0);
  qc2qr_autoconv_tend.set(nc2nr_autoconv_tend == 0 && context, 0);

}

} // namespace p3
} // namespace scream

#endif
