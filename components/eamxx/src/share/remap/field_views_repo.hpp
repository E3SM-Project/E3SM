#ifndef EAMXX_FIELD_VIEWS_REPO_HPP
#define EAMXX_FIELD_VIEWS_REPO_HPP

#include "share/field/field.hpp"
#include "share/core/eamxx_types.hpp"

#include <map>
#include <vector>

namespace scream
{

/*
 * A repository of device views for a list of fields
 *
 * Remappers that process all registered fields in a single kernel need a way
 * to access the views of all fields from inside the kernel. This class stores
 * the views (up to rank 5) of a list of fields in a device view of containers,
 * so that the i-th entry of 'views'/'cviews' can be accessed in a kernel.
 *
 * Fields that are subfields of other fields may have a subview info that
 * changes over time (see Field::subfield). For those, the views must be
 * re-extracted before every use, which is what update_subfields does.
 */

class FieldViewsRepo
{
public:
  using device_type = Field::device_t;
  using KT          = KokkosTypes<device_type>;

  template<typename T,int N>
  using view_Nd = typename KT::template view_ND<T,N>;
  template<typename T>
  using view_1d = view_Nd<T,1>;

  template<typename T>
  struct ViewsContainer {
    view_Nd<T,1>   v1d;
    view_Nd<T,2>   v2d;
    view_Nd<T,3>   v3d;
    view_Nd<T,4>   v4d;
    view_Nd<T,5>   v5d;
  };

  using views_t   = view_1d<ViewsContainer<Real>>;
  using cviews_t  = view_1d<ViewsContainer<const Real>>;
  using hviews_t  = typename views_t::host_mirror_type;
  using hcviews_t = typename cviews_t::host_mirror_type;

  // Store the views of the input fields. Only the const views are stored for read-only fields.
  // Fields that are subfields of other fields (but not of subfields) are tracked,
  // so that update_subfields can detect changes in their subview info.
  void init (const std::vector<Field>& fields);

  // Check if the subview info of any tracked subfield changed since init (or the last update)
  bool subfields_info_have_changed (const std::vector<Field>& fields) const;

  // Re-extract the views of the subfields whose subview info changed
  void update_subfields (const std::vector<Field>& fields);

  views_t   views;
  cviews_t  cviews;

  hviews_t  h_views;
  hcviews_t h_cviews;

private:

  void set_views (const int i, const Field& f);

  // For each subfield, we store the last value of subview info, to check
  // if anything changed (only matters for 'dynamic' subfields).
  // Note: here "dynamic" is in the sense explained in Field::subfield
  std::map<int,SubviewInfo> m_subfield_info;
};

// Reinterpret a view of Real as a view of NewValueT (e.g., a pack of Real).
// NOTE: one could deduce from OldViewT whether the return type
//       should have a const value type. But the code is a bit tedious,
//       so we'll just force to call this as pack_view<const T>(v).
template<typename NewValueT, typename OldViewT>
KOKKOS_INLINE_FUNCTION
FieldViewsRepo::view_Nd<NewValueT,OldViewT::rank>
pack_view (const OldViewT& v)
{
  constexpr int N = OldViewT::rank;
  Kokkos::LayoutRight kl;
  for (int i=0; i<N; ++i) {
    kl.dimension[i] = v.extent(i);
  }
  using OldValueT = typename OldViewT::traits::value_type;
  auto pack_size = sizeof(NewValueT) / sizeof(OldValueT);
  kl.dimension[N-1] /= pack_size;
  auto ptr = reinterpret_cast<NewValueT*>(v.data());

  // Naively, you would return tmp directly. However, if v is a
  // subview of another view along the second extent, you need
  // to fix the stride, since the ctor from ptr+layout would have
  // no knowledge about that.
  FieldViewsRepo::view_Nd<NewValueT,N> tmp (ptr,kl);
  auto vm = tmp.impl_map();
  vm.m_impl_offset.m_stride = v.impl_map().stride_0() / pack_size;
  return FieldViewsRepo::view_Nd<NewValueT,N> (tmp.impl_track(),vm);
}

} // namespace scream

#endif // EAMXX_FIELD_VIEWS_REPO_HPP
