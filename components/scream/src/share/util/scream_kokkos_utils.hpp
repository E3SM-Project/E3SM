#ifndef SCREAM_KOKKOS_UTILS_HPP
#define SCREAM_KOKKOS_UTILS_HPP

#include <share/scream_kokkos_meta.hpp>
#include <share/util/scream_std_meta.hpp>

#include <Kokkos_Core.hpp>

#include <cassert>
#include <type_traits>

// This file should not be merged with share/scream_kokkos_meta.hpp.
// That file contains functionalities that *ideally* should be in kokkos
// itself, and we can foresee being in kokkos in the future. Hence, that
// file may be gone at some point (in fact, it may be gone by the time
// you read this comment). This file, instead, contains functionalities
// that are probably not generic enough to appear in kokkos any time soon
// (or ever), and are more app-specific.

namespace scream {
namespace util {

// Retrieve the underlying scalar type from a MD array type
template<typename MDArrayType>
struct ScalarType {
private:
  using no_ref = typename std::remove_reference<MDArrayType>::type;
  using no_arr = typename std::remove_all_extents<no_ref>::type;
public:
  using type = typename remove_all_pointers<no_arr>::type;
};

template<typename DataType>
struct GetRanks {
private:
  using dimension = typename Kokkos::View<DataType>::traits::dimension;
public:
  enum : int { rank = dimension::rank };
  enum : int { rank_dynamic = dimension::rank_dynamic };
};

template<typename DataTypeOut, typename DataTypeIn, typename... Props>
typename std::enable_if<GetRanks<DataTypeOut>::rank_dynamic==0,
                        typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>>>::type
reshape (Kokkos::View<DataTypeIn,Props...> view_in) {
  typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>> view_out(view_in.data());
  assert (view_in.size()==view_out.size());
  return view_out;
}

template<typename DataTypeOut, typename DataTypeIn, typename... Props>
typename std::enable_if<GetRanks<DataTypeOut>::rank_dynamic==1,
                        typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>>>::type
reshape (Kokkos::View<DataTypeIn,Props...> view_in,
         const int dim0) {
  typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>> view_out(view_in.data(),dim0);
  assert (view_in.size()==view_out.size());
  return view_out;
}

template<typename DataTypeOut, typename DataTypeIn, typename... Props>
typename std::enable_if<GetRanks<DataTypeOut>::rank_dynamic==2,
                        typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>>>::type
reshape (Kokkos::View<DataTypeIn,Props...> view_in,
         const int dim0, const int dim1) {
  typename ko::Unmanaged<Kokkos::View<DataTypeOut,Props...>> view_out(view_in.data(),dim0,dim1);
  assert (view_in.size()==view_out.size());
  return view_out;
}

} // namespace util
} // namespace scream

#endif // SCREAM_KOKKOS_UTILS_HPP
