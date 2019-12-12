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

template <typename PackType, typename ExpType>
KOKKOS_INLINE_FUNCTION
PackType bfb_pow (const PackType& a, const ExpType e) {
  PackType b;
  for (int s = 0; s < PackTraits<PackType>::pack_length; ++s) {
    b[s] = std::pow(a[s], e);
  }

#ifdef CUDA_BUILD
  zero_ulp_n(b, 25, a);
#endif
  return b;
}

} // namespace Homme

#endif // HOMME_BFB_UTILS_HPP
