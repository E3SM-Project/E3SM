#include "physics_share.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

#include <random>

namespace scream {

namespace physics {

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
