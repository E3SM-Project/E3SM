#ifndef EKAT_MD_ARRAY_HPP
#define EKAT_MD_ARRAY_HPP

#include <array>

/*
 * A template alias and utilities for potentially nested std::arrays
 * 
 * When using compile-time arrays, std::array should be preferred over
 * C-style arrays. However, syntax become heavy when we have nested
 * std::arrays. The template alias defined here allows one to write
 * 
 *    md_array<double,2,3,4> a;
 * 
 * instead of
 * 
 *    std::array<std::array<std::array<double,2>,3>,4> a;
 * 
 * Furthermore, we provide a free function 'data', that extract the
 * pointer to the actual data. Notice that, if a is defined as above
 * calling 'a.data()' would not return a pointer to double, but
 * a pointer to std::array<std::array<double,2>,3>, which in some cases
 * is not what we want. Instead, calling `data(a)` returns a double*,
 * which is what you'd get with &a[0][0][0], but having to know
 * the actual rank of the md array.
 */

namespace ekat {

namespace util {

// Helper for MD std::array
template<typename T, std::size_t M, std::size_t...N>
struct md_array_helper {
  using type = std::array<typename md_array_helper<T,N...>::type,M>;
};

template<typename T, std::size_t M>
struct md_array_helper<T,M> {
  using type = std::array<T,M>;
};

// Shorter name for an md type
template<typename T, std::size_t... N>
using md_array = typename md_array_helper<T,N...>::type;

// Meta-utility helper structure
template<typename T>
struct get_md_data;

// Recursion case, where std::array's value type is itself an std::array
template<typename T, std::size_t M, std::size_t N>
struct get_md_data<std::array<std::array<T,N>,M>> {
  using type = std::array<std::array<T,N>,M>;
  using value_type = typename get_md_data<std::array<T,N>>::value_type;

  static value_type* data (type& a) {
    return get_md_data<std::array<T,N>>::data(a[0]);
  }
  static const value_type* data (const type& a) {
    return get_md_data<std::array<T,N>>::data(a[0]);
  }
};

// Base case, where std::array's value type is not itself an std::array
template<typename T, std::size_t N>
struct get_md_data<std::array<T,N>> {
  using type = std::array<T,N>;
  using value_type = T;

  static T* data (type& a) {
    return a.data();
  }

  static const T* data (const type& a) {
    return a.data();
  }
};

// Free function utilities to extract a pointer to the flattened md array data
template<typename T, std::size_t N>
typename get_md_data<std::array<T,N>>::value_type*
data (std::array<T,N>& a) {
  return get_md_data<std::array<T,N>>::data(a);
}

template<typename T, std::size_t N>
const typename get_md_data<std::array<T,N>>::value_type*
data (const std::array<T,N>& a) {
  return get_md_data<std::array<T,N>>::data(a);
}

} // namespace util

} // namespace ekat

#endif // EKAT_MD_ARRAY_HPP
