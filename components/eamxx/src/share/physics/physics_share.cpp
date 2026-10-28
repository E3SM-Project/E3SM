#include "physics_share.hpp"
#include "share/core/eamxx_config.hpp"
#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"

#include <ekat_assert.hpp>

#include <cmath>
#include <random>

namespace scream {
namespace physics {

SolarGeometry solar_geometry(const util::TimeStamp& timestamp, const ekat::ParameterList& params,
                             int orbital_year_setting, double eccentricity,
                             double obliquity, double perihelion) {
  const auto calday = timestamp.frac_of_year_in_days() + 1;
  const auto solar_ymd = timestamp.get_year() * 10000 + timestamp.get_month() * 100 + timestamp.get_day();
  const auto solar_tod = (calday - std::floor(calday)) * 86400;
  const bool has_coupled_geometry = shr_orb_has_geometry_callback_c2f();
  const bool coupled_eot = shr_orb_use_equation_of_time_c2f();
  const bool has_local_eot = params.isParameter("use_equation_of_time");
  const bool local_eot = has_local_eot ? params.get<bool>("use_equation_of_time") : true;
  EKAT_REQUIRE_MSG(!has_coupled_geometry || !has_local_eot || local_eot == coupled_eot,
                  "Solar equation-of-time setting conflicts with the coupled driver.");
  const bool use_eot = has_coupled_geometry ? coupled_eot : local_eot;
  double delta, eccf, eot_seconds = 0;
  const bool coupled_orbit = has_coupled_geometry
    ? shr_orb_coupled_geometry_c2f(solar_ymd, solar_tod, &delta, &eccf, &eot_seconds)
    : use_leap_year() && shr_orb_coupled_decl_c2f(solar_ymd, solar_tod, &delta, &eccf);
  EKAT_REQUIRE_MSG(!has_coupled_geometry || coupled_orbit, "Coupled solar geometry is unavailable.");
  if (!coupled_orbit) {
    double obliqr, lambm0, mvelpp;
    auto orbital_year = orbital_year_setting;
    auto eccen = eccentricity;
    auto obliq = obliquity;
    auto mvelp = perihelion;
    if (eccen >= 0 && obliq >= 0 && mvelp >= 0) {
      orbital_year = shr_orb_undef_int_c2f;
    } else if (orbital_year < 0) {
      orbital_year = timestamp.get_year();
    }
    if (use_leap_year() && orbital_year_setting < 0 && orbital_year != shr_orb_undef_int_c2f) {
      const auto year_fraction = (calday - 1) / timestamp.days_in_curr_year();
      shr_orb_params_continuous_c2f(orbital_year, year_fraction, &eccen, &obliq,
                                  &mvelp, &obliqr, &lambm0, &mvelpp);
    } else {
      shr_orb_params_c2f(&orbital_year, &eccen, &obliq, &mvelp, &obliqr, &lambm0, &mvelpp);
    }
    if (use_leap_year() && orbital_year_setting >= 0 && orbital_year != shr_orb_undef_int_c2f) {
      shr_orb_decl_fixed_year_c2f(calday, timestamp.days_in_curr_year(),
                                util::is_leap_year(timestamp.get_year()-1) ? 366 : 365,
                                eccen, mvelpp, lambm0, obliqr, &delta, &eccf, &eot_seconds, use_eot);
    } else if (use_leap_year()) {
      const util::TimeStamp epoch(2000, 1, 1, 0, 0, 0);
      const util::TimeStamp midnight(timestamp.get_year(), timestamp.get_month(), timestamp.get_day(), 0, 0, 0);
      const auto elapsed_days = midnight.days_from(epoch) + (calday - std::floor(calday));
      if (use_eot) {
        shr_orb_decl_gregorian_eot_c2f(elapsed_days, eccen, mvelpp, lambm0,
                                     obliqr, &delta, &eccf, &eot_seconds);
      } else {
        shr_orb_decl_gregorian_c2f(elapsed_days, eccen, mvelpp, lambm0, obliqr, &delta, &eccf);
      }
    } else if (use_eot) {
      shr_orb_decl_eot_c2f(calday, eccen, mvelpp, lambm0, obliqr, &delta, &eccf, &eot_seconds);
    } else {
      shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0, obliqr, &delta, &eccf);
    }
  }
  return SolarGeometry(calday, delta, eccf, eot_seconds, use_eot);
}

double SolarGeometry::cosine_zenith(double latitude_radians, double longitude_radians,
                                  double averaging_seconds) const {
  EKAT_REQUIRE_MSG(std::isfinite(averaging_seconds) && averaging_seconds >= 0,
                  "Solar cosine averaging interval must be finite and nonnegative.");
  return m_use_equation_of_time
    ? shr_orb_cosz_eot_c2f(m_calendar_day, latitude_radians, longitude_radians, m_declination,
                         averaging_seconds, m_equation_of_time)
    : shr_orb_cosz_c2f(m_calendar_day, latitude_radians, longitude_radians, m_declination, averaging_seconds);
}

//
// A C++ interface to physics-related fortran calls and vice versa
//

// Cuda implementations of std math routines are not necessarily BFB
// with the host.
template <typename ScalarT, typename DeviceT>
struct CudaWrap
{
  using Scalar = ScalarT;
  using RangePolicy = typename ekat::KokkosTypes<DeviceT>::RangePolicy;

  static Scalar pow(Scalar base, Scalar exp)
  {
    Scalar result;
    RangePolicy policy(0,1);
    Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const Int&, Scalar& value) {
        value = std::pow(base, exp);
    }, result);

    return result;
  }

#define cuda_wrap_single_arg(wrap_name, func_call)                            \
static Scalar wrap_name(Scalar input) {                                       \
  Scalar result;                                                              \
  RangePolicy policy(0,1);                                                    \
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const Int&, Scalar& value) {  \
    value = func_call(input);                                                 \
  }, result);                                                                 \
  return result;                                                              \
}

  cuda_wrap_single_arg(gamma, std::tgamma)
  cuda_wrap_single_arg(sqrt, std::sqrt)
  cuda_wrap_single_arg(cbrt, std::cbrt)
  cuda_wrap_single_arg(log, std::log)
  cuda_wrap_single_arg(log10, std::log10)
  cuda_wrap_single_arg(exp, std::exp)
  cuda_wrap_single_arg(expm1, std::expm1)
  cuda_wrap_single_arg(tanh, std::tanh)
  cuda_wrap_single_arg(cos, std::cos)
  cuda_wrap_single_arg(erf, std::erf)

#undef cuda_wrap_single_arg
};

extern "C" {

Real scream_pow(Real base, Real exp)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::pow(base, exp);
#else
  return std::pow(base, exp);
#endif
}

Real scream_gamma(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::gamma(input);
#else
  return std::tgamma(input);
#endif
}

Real scream_cbrt(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::cbrt(input);
#else
  return std::cbrt(input);
#endif
}

Real scream_sqrt(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::sqrt(input);
#else
  return std::sqrt(input);
#endif
}

Real scream_log(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::log(input);
#else
  return std::log(input);
#endif
}

Real scream_log10(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::log10(input);
#else
  return std::log10(input);
#endif
}

Real scream_exp(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::exp(input);
#else
  return std::exp(input);
#endif
}

Real scream_expm1(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::expm1(input);
#else
  return std::expm1(input);
#endif
}

Real scream_tanh(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::tanh(input);
#else
  return std::tanh(input);
#endif
}

Real scream_cos(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::cos(input);
#else
  return std::cos(input);
#endif
}

Real scream_erf(Real input)
{
#ifdef EAMXX_ENABLE_GPU
  return CudaWrap<Real, DefaultDevice>::erf(input);
#else
  return std::erf(input);
#endif
}

} // end extern "C"

} // namespace physics
} // namespace scream
