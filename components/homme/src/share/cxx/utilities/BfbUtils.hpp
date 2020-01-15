#ifndef HOMME_BFB_UTILS_HPP
#define HOMME_BFB_UTILS_HPP

#include "PackTraits.hpp"
#include <Kokkos_Core.hpp>

namespace Homme
{
// Set the least significant n bits of the fraction to 0. Use this for
// BFB-testing between GPU and host. Zero the n least significant bits of y in
//     y = pow(x, e),
// for example.
//   Table 8 of Appendix E of the Cuda programming guide gives the accuracy of
// various functions. However, we also have to account for carries in the
// calculations. Thus, even if an answer is accurate to 2 ULP, carries can make
// more bits differ. Thus, zeroing (or any sort of rounding) is not assured to
// work. However, you can reduce the probability of a diff to a small number by
// choosing n to be large. Full double precision range is 53 bits; choose
// something like n = 20 to get a carry-based diff every 1 in a million calls,
// still leaving 33 bits.
//   To remain sensitive to the full bits in the argument to the non-BFB
// function, such as x in pow above, pass replace = x. The n least significant
// bits of the fraction part of x is used to fill the 0s introduced in this
// function. Thus, the final result z of
//     y = pow(x, e);
//     z = zeroulpn(y, n, x)
// is sensitive to all bits of x.
KOKKOS_INLINE_FUNCTION
double zeroulpn (double a, const int n, double replace) {
  const auto f = static_cast<double>(1LL << 54);
  int e;
  a = std::frexp(a, &e);
  const int sign = a >= 0 ? 1 : -1;
  a = std::abs(a);
  long long int x = f*a;
  x >>= n;
  x <<= n;
  { // Replace the n least significant bits with those from the fraction part of
    // 'replace'.
    const long long int ones = 0xffffffffffffffff;
    auto mask = ones;
    mask >>= n;
    mask <<= n;
    mask = ones & (mask ^ ones);
    int er;
    replace = std::abs(std::frexp(replace, &er));
    long long int r = f*replace;
    x = x | (mask & r);
  }
  a = x/f;
  return std::ldexp(sign*a, e);
}

template <typename PackType>
KOKKOS_INLINE_FUNCTION
void zero_ulp_n (PackType& a, const int n, const PackType& replace) {
  for (int s = 0; s < PackTraits<PackType>::pack_length; ++s) {
    a[s] = zeroulpn(a[s], n, replace[s]);
  }
}

template<typename ScalarType>
KOKKOS_INLINE_FUNCTION
ScalarType int_pow (const ScalarType val, int k) {
  if (k<=0) {
    Kokkos::abort("int_pow implemented only for the case k>=1.\n");
  }
  if (k==1) {
    return val;
  } else if (k%2 == 0) {
    return int_pow (val*val, k/2);
  } else {
    return val*int_pow(val,k-1);
  }
}

template <typename ScalarType, typename ExpType>
KOKKOS_INLINE_FUNCTION
ScalarType bfb_pow_impl (ScalarType val, const ExpType e) {
#ifdef CUDA_BUILD
  if (val==ScalarType(0)) {
    if (e==ExpType(0)) {
      Kokkos::abort("Cannot do 0^0, sorry man.\n");
    }
    return ScalarType(0);
  } else if (val>=ScalarType(1.5)) {
    int k=0;
    constexpr ScalarType two(2.0);
    constexpr ScalarType sixteen(16.0);
    while (val>=sixteen){
      k += 4;
      val /= sixteen;
    }

    while (val>=ScalarType(1.5)){
      ++k;
      val /= two;
    }

    // 2^e approx with 4th order Taylor expansion
    ScalarType two_e = 1 + e*(1 + (e-1)/2*(1 + (e-2)/3*(1 + (e-3)/4)));
    
    return int_pow(two_e,k)*bfb_pow_impl(val,e);
  } else if (val<=ScalarType(0.5)) {
    int k=0;
    constexpr ScalarType two(2.0);
    constexpr ScalarType sixteen(16.0);
    constexpr ScalarType half(0.5);
    constexpr ScalarType sixteenth(1.0/16.0);
    while (val<=sixteenth){
      k += 4;
      val *= sixteen;
    }

    while (val<=half){
      ++k;
      val *= two;
    }

    // 2^e approx with 4th order Taylor expansion
    ScalarType two_e = 1 + e*(1 + (e-1)/2*(1 + (e-2)/3*(1 + (e-3)/4)));
    
    return bfb_pow_impl(val,e)/int_pow(two_e,k);
  } else {
    constexpr int nmax = 10;

    ScalarType x0 = val - ScalarType(1.0);
    ScalarType a0 = e;
    ScalarType y (0.0);

    ScalarType x(1.0);
    ScalarType a(1.0);
    for (int n=0; n<=nmax; ++n) {
      y += a*x;
      x *= x0;
      a *= (a0-n)/(n+1);
    }
    return y;
  }
#else
  return std::pow(val,e);
#endif
}

template <typename PackType, typename ExpType>
KOKKOS_INLINE_FUNCTION
PackType bfb_pow (const PackType& a, const ExpType e) {
  PackType b;
  for (int s = 0; s < PackTraits<PackType>::pack_length; ++s) {
    b[s] = bfb_pow_impl(a[s], e);
  }
  return b;
}

template <>
KOKKOS_INLINE_FUNCTION
double bfb_pow<double,double> (const double& a, const double e) {
  return bfb_pow_impl(a,e);
}

} // namespace Homme

#endif // HOMME_BFB_UTILS_HPP
