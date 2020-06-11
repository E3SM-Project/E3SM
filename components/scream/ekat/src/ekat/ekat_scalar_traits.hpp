#ifndef EKAT_SCALAR_TRAITS_HPP
#define EKAT_SCALAR_TRAITS_HPP

#include "scream_kokkos.hpp"

#include <limits>
#include <type_traits>

namespace scream {

template<typename T>
struct ScalarTraits {

  static_assert (std::is_floating_point<T>::value, "Error! Template parameter 'T' in ScalarTraits must be a numeric type.\n");

  static const T quiet_NaN () {
    static bool inited = false;
    static T nan;
    if (!inited) {
      nan = get_quiet_NaN();
      inited = true;
    }
    return nan;
  }

private:

  static const T get_quiet_NaN () {
#ifdef KOKKOS_ENABLE_CUDA
    auto policy = Kokkos::RangePolicy<ExeSpaceType>(0,1);
    Kokkos::View<double[1],Kokkos::Cuda> d_nan("");
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int) {
      d_nan(0) = nan("");
    });
    auto h_nan = Kokkos::create_mirror_view(d_nan);
    Kokkos::deep_copy(h_nan,d_nan);
    return h_nan(0);
#else
    return std::numeric_limits<T>::quiet_NaN();
#endif
  }
};


} // namespace scream

#endif // EKAT_SCALAR_TRAITS_HPP
