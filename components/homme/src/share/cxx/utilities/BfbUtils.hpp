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
    const unsigned long long int ones = 0xffffffffffffffff;
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
ScalarType int_pow (ScalarType val, int k) {
  constexpr int max_shift = 30;
  if (k<0) {
    printf ("k = %d\n",k);
    Kokkos::abort("int_pow implemented only for k>=0.\n");
  }

  ScalarType y(1);
  if (k==0) {
    return y;
  }

  for (int shift=0; shift<max_shift; ++shift) {
    int pow2 = 1<<shift;
    if (k & pow2) {
      y *= val;
    }
    val *= val;
  }
  return y;
}

template <typename ScalarType>
KOKKOS_INLINE_FUNCTION
ScalarType power_of_two (const ScalarType e) {
  // Note: this function is tailored (or taylored...eheh), for -1<e<1.5

  constexpr ScalarType l2 = 0.693147180559945;

  const ScalarType l2e = l2*e;

  ScalarType y (1);

  // Note: if order = 0, y=1, which is the initial value.
  constexpr int order = 10;
  for (int n=order; n>0; --n) {
    y *= l2e/n;
    y += 1.0;
  }
  return y;
}

// Note: force type to not be reference, so we can
//       recycle inputs during calculations.
template <typename ScalarType, typename ExpType>
KOKKOS_INLINE_FUNCTION
typename
std::enable_if<!std::is_reference<ScalarType>::value &&
               !std::is_reference<ExpType>::value, ScalarType>::type
bfb_pow_impl (ScalarType val, ExpType e) {
#ifdef HOMMEXX_ENABLE_GPU
  // Note: this function is tailored (or taylored...eheh)
  // for -1<e<1.5 and 0.001 < b < 1e6
  if (val<ScalarType(0)) {
    Kokkos::abort("Cannot take powers of negative numbers.\n");
  }
  if (e<-1.0 || e>1.5) {
    Kokkos::abort("bfb_pow x^a impl-ed with only -1.0<a<1.5 in mind.\n");
  }
  if (val==ScalarType(0)) {
    if (e==ExpType(0)) {
      Kokkos::abort("Cannot do 0^0, sorry man.\n");
    }
    return ScalarType(0);
  }

  if (e<0) {
    e = -e;
    val = 1/val;
  }

  ScalarType factor = 1.0;
  if (val>=ScalarType(1.5)) {
    int k=0;
    while (val>=16.0){
      k += 4;
      val /= 16.0;
    }

    while (val>=1.5){
      ++k;
      val /= 2.0;
    }

    factor = int_pow(power_of_two(e),k);
  } else if (val<=0.5) {
    int k=0;
    while (val<=(1.0/16)){
      k += 4;
      val *= 16.0;
    }

    while (val<=0.5){
      ++k;
      val *= 2.0;
    }

    factor = 1.0/int_pow(power_of_two(e),k);
  }

  ScalarType x = val - ScalarType(1.0);
  ScalarType y (1.0);

  // (1+x)^a = 1+ax+a*(a-1)/2 x^2 + a*(a-1)*(a-2)/6 x^3 +...
  //         = 1+ax*(1+(a-1)/2 x *( 1+(a-2)/3 x * (1+...)))
  // Note: if order = 0, y=1, which is the initial value.
  constexpr int order = 5;
  for (int n=order; n>=1; --n) {
    y *= ((e-(n-1))/n) * x;
    y += 1.0;
  }

  return y*factor;
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
