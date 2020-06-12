#ifndef EKAT_SCALAR_TRAITS_HPP
#define EKAT_SCALAR_TRAITS_HPP

#include "scream_kokkos.hpp"

#include <limits>
#include <climits>
#include <type_traits>

#ifdef KOKKOS_ENABLE_CUDA
#include <math_constants.h>
#endif

namespace scream {

template<typename T>
struct ScalarTraits {

  static_assert (std::is_arithmetic<T>::value, "Error! Template parameter 'T' in ScalarTraits must be a numeric type.\n");

  KOKKOS_INLINE_FUNCTION
  static const T quiet_NaN () {
    scream_kassert_msg(std::is_floating_point<T>::value, "Error! NaN is only available for floating point types.\n");
#ifdef __CUDA_ARCH__
    if (std::is_same<T,float>::value) {
      return CUDART_NAN_F;
    } else if (std::is_same<T,float>::value) {
      return CUDART_NAN;
    } else {
      scream_kerror_msg ("Error! No NaN provided for this floating point type.\n");
    }
#else
    return std::numeric_limits<T>::quiet_NaN();
#endif
  }

  KOKKOS_INLINE_FUNCTION
  static const T invalid () {
    // For a floating point, return NaN. For an integer, return largest possible number
    if (std::is_floating_point<T>::value) {
      return quiet_NaN();
    } else {
      // If cuda supported numeric_limits, we would not need this
      if (std::is_same<T,int>::value) {
        return INT_MAX;
      } else if (std::is_same<T,unsigned int>::value) {
        return UINT_MAX;
      } else if (std::is_same<T,long>::value) {
        return LONG_MAX;
      } else if (std::is_same<T,unsigned long>::value) {
        return ULONG_MAX;
      } else if (std::is_same<T,long long>::value) {
        return LLONG_MAX;
      } else if (std::is_same<T,unsigned long long>::value) {
        return ULLONG_MAX;
      }
    }
  }
};


} // namespace scream

#endif // EKAT_SCALAR_TRAITS_HPP
