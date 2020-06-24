#ifndef EKAT_SCALAR_TRAITS_HPP
#define EKAT_SCALAR_TRAITS_HPP

#include "scream_kokkos.hpp"

#include <limits>
#include <climits>
#include <type_traits>
#include <typeinfo>

#ifdef KOKKOS_ENABLE_CUDA
#include <math_constants.h>
#endif

namespace scream {

template<typename T>
struct ScalarTraits {

  using raw_type = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

  static_assert (std::is_arithmetic<raw_type>::value, "Error! Template parameter 'T' in ScalarTraits must be a numeric type.\n");

  KOKKOS_INLINE_FUNCTION
  static const raw_type quiet_NaN () {
    scream_kassert_msg(std::is_floating_point<raw_type>::value,
                       "Error! NaN is only available for floating point types.\n");
#ifdef __CUDA_ARCH__
    if (std::is_same<raw_type,float>::value) {
      return CUDART_NAN_F;
    } else if (std::is_same<raw_type,double>::value) {
      return CUDART_NAN;
    } else {
      scream_kerror_msg ("Error! No NaN provided for this floating point type.\n");
      // Silence compiler warning
      return 0;
    }
#else
    return std::numeric_limits<raw_type>::quiet_NaN();
#endif
  }

  KOKKOS_INLINE_FUNCTION
  static const raw_type invalid () {
    // For a floating point, return NaN. For an integer, return largest possible number
    raw_type val(0);

    if (std::is_floating_point<raw_type>::value) {
      val = ScalarTraits<raw_type>::quiet_NaN();
    } else {
      // If cuda supported numeric_limits, we would not need all these ifs
      if (std::is_same<raw_type,int>::value) {
        val = static_cast<raw_type>(INT_MAX);
      } else if (std::is_same<raw_type,unsigned int>::value) {
        val = static_cast<raw_type>(UINT_MAX);
      } else if (std::is_same<raw_type,long>::value) {
        val = static_cast<raw_type>(LONG_MAX);
      } else if (std::is_same<T,unsigned long>::value) {
        val = static_cast<raw_type>(ULONG_MAX);
      } else if (std::is_same<raw_type,long long>::value) {
        val = static_cast<raw_type>(LLONG_MAX);
      } else if (std::is_same<raw_type,unsigned long long>::value) {
        val = static_cast<raw_type>(ULLONG_MAX);
      }
    }
    return val;
  }
};


} // namespace scream

#endif // EKAT_SCALAR_TRAITS_HPP
