#ifndef EKAT_MD_ARRAY_HPP
#define EKAT_MD_ARRAY_HPP

#include <array>

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

// Meta-utility to get inner-most value type
template<typename T>
struct get_md_data;

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
