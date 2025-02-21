#ifndef SCREAM_FIELD_ALLOC_PROP_HPP
#define SCREAM_FIELD_ALLOC_PROP_HPP

#include "share/field/field_layout.hpp"
#include "share/eamxx_types.hpp"

#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/ekat_assert.hpp"

#include <vector>
#include <memory>

namespace scream
{

/*
 *  Small structure holding a few properties of the field allocation
 *
 *  This class allows to keep track of a field allocation properties.
 *  There are three needs for these info: allocation, reshaping, and subviews.
 *    - At allocation time, we need to know the alloc size, but this might
 *      be *larger* than what the field layout would suggest, in case
 *      padding was used to accommodate use of the fields with packs.
 *
 *      For instance, say you declare a field with dimensions (3,10).
 *      In the field manager, this would be stored as a 1d field: View<Real*>.
 *      Let's say you have one class that uses the field as View<Real**>,
 *      and another class that used the field in a packed way, that is,
 *      as View<Pack<Real,4>**>. The second view will have dimensions (3,3),
        in terms of its value type (i.e., Pack<Real,4>).
 *      This means the field needs an allocation bigger than 30 Real's, namely
 *      3*4*3=36 Real's. This class does the book-keeping for the allocation size,
 *      so that 1) the field can be allocated with enough memory to accommodate
 *      all requests, 2) customers of the field can check what the allocation
 *      is, so that they know whether there is padding in the field, and
 *      3) query whether the allocation is compatible with a given value type.
 *      The last point allows one to check if the view can be reshaped
 *      as, say, View<Pack<Real,16>> or not. Notice that a field can be
 *      reshaped as, say, View<Pack<Real,16>> even if no request for
 *      such large pack was made. E.g., if the field dimensions are
 *      (10,128), then the last dim allows packing with pack size 16.
 *    - Fields are stored as 1d views (basically, just a pointer), so
 *      users need to call the get_view<T>() to use it with
 *      a different data type. E.g., to use the field as a 2d field,
 *      with Pack<Real,8> as scalar type, one needs to call
 *        auto v = f.get_view<Pack<Real,8>**>()
 *      The alloc props are queried to establish 1) whether the allocation
 *      is compatible with a pack-8 scalar type, and 2) how the dimensions
 *      of the resulting view should be.
 *    - It is possible (with some limitations) to have a Field being a
 *      "subview" of a bigger field. E.g., tracers are allocated as a big
 *      array Q, but one may just be interested in a single one (e.g., qv),
 *      yet still use it as a Field. This requires some care during the
 *      get_view method, so the field alloc props can be used
 *      to detect whether we are in this scenario.
 *
 *  Note: at every request for a new value_type, this class checks the
 *        underlying scalar_type. We ASSUME that the dimensions in the
 *        field layout refer to that scalar_type.
 */

// Helper struct to store info related to the subviewing process
struct SubviewInfo {
  SubviewInfo() = default;
  SubviewInfo(const int dim, const int slice, const int extent,
              const bool is_dynamic)
      : dim_idx(dim), slice_idx(slice), dim_extent(extent),
        dynamic(is_dynamic) {}
  // multi-slice subview across contiguous indices
  SubviewInfo(const int dim, const int slice_beg, const int slice_end,
              const int extent)
      : dim_idx(dim), slice_idx(slice_beg), slice_idx_end(slice_end),
        dim_extent(extent) {}
  SubviewInfo(const SubviewInfo&) = default;
  SubviewInfo& operator=(const SubviewInfo&) = default;

  // Dimension along which slicing happened
  int dim_idx = -1;
  // in the case of single slice, this is the slice index
  // in the case of multi-slice, this is the beginning index
  int slice_idx = -1;
  // ending slice index for multi-slice (remains == -1 if single-slice subview)
  int slice_idx_end = -1;
  // Extent of dimension $dim_idx
  int dim_extent = -1;
  // Whether this is a dynamic subview (slice_idx can change)
  bool dynamic = false;
};

inline bool operator== (const SubviewInfo& lhs, const SubviewInfo& rhs) {
  return lhs.dim_idx==rhs.dim_idx &&
         lhs.slice_idx==rhs.slice_idx &&
         //  slice_idx_end == -1 for single slice, and the ending index when
         // it's a multi-slice
         lhs.slice_idx_end==rhs.slice_idx_end &&
         lhs.dim_extent==rhs.dim_extent &&
         lhs.dynamic==rhs.dynamic;
}

class FieldAllocProp {
public:
  using layout_type      = FieldLayout;
  using layout_ptr_type  = std::shared_ptr<const layout_type>;

  FieldAllocProp (const int scalar_size);
  FieldAllocProp (const FieldAllocProp&) = default;

  FieldAllocProp& operator= (const FieldAllocProp&);

  // Return allocation props of an unmanaged subview of this field
  // at entry k along dimension idim. If dynamic = true, then it will be
  // possible to change the entry index (k) at runtime.
  FieldAllocProp subview (const int idim, const int k, const bool dynamic) const;

  // multi-slice subview over contiguous indices
  FieldAllocProp subview (const int idim, const int k_beg, const int k_end) const;

  // Request allocation able to accommodate a pack of ScalarType of the given pack size
  void request_allocation (const int pack_size = 1);

  // Request allocation able to accommodate all the alloc props in src.
  // Note: src does not need to be committed yet.
  void request_allocation (const FieldAllocProp& src);

  // Locks the allocation properties, preventing furter value types requests
  void commit (const layout_type& layout);

  // For dynamic subfield, reset the slice index
  void reset_subview_idx (const int idx);

  // ---- Getters ---- //

  // Whether commit has been called
  bool is_committed () const { return m_committed; }

  // Whether this field is a subfield of another
  bool is_subfield () const { return m_subview_info.dim_idx>=0; }

  // Whether this field is a dynamic subfield of another
  bool is_dynamic_subfield () const { return m_subview_info.dynamic; }

  // Get the overall allocation size (in Bytes)
  long long get_alloc_size () const;

  // Get number of m_scalar_type_size-sized scalars in the allocation.
  long long get_num_scalars () const { return get_alloc_size () / m_scalar_type_size; }

  // Get the largest pack size that this allocation was requested to accommodate
  int get_largest_pack_size () const;

  // Wether this allocation is contiguous
  bool contiguous () const { return m_contiguous; }

  // Size of the last extent in the alloction (i.e., number of scalars in it)
  int  get_last_extent () const;
  int  get_padding () const;

  // Return the slice information of this subview allocation.
  const SubviewInfo& get_subview_info () const { return m_subview_info; }

  // Whether the allocation properties are compatible with ValueType
  template<typename ValueType>
  bool is_compatible () const;

protected:

  // The FieldLayout associated to this allocation
  layout_type     m_layout;

  // The list of requested value types for this allocation
  std::vector<int>    m_value_type_sizes;

  // The size of the scalar type
  int   m_scalar_type_size;

  // The largest pack size that was requested
  int   m_pack_size_max;

  // The extent along the last dimension for this allocation
  int   m_last_extent;

  // The total size of this allocation.
  long long   m_alloc_size;

  // If this allocation is a subview of another, store info about the subview
  // process (see SubviewInfo documentation for details)
  SubviewInfo         m_subview_info;

  // Whether the allocation is contiguous
  bool  m_contiguous;

  // Whether commit was called (i.e., no more value type requests allowed)
  bool  m_committed;
};

// ================================= IMPLEMENTATION ================================== //

template<typename ValueType>
bool FieldAllocProp::is_compatible () const {
  using ekat::ScalarTraits;

  constexpr int sts = sizeof(typename ScalarTraits<ValueType>::scalar_type);
  constexpr int vts = sizeof(ValueType);
  constexpr int vlen = vts / sts;

  // Allocation is compatible with ValueType if
  //   - the scalar type of ValueType has the same size as m_scalar_type_size
  //   - the size of ValueType must divide the last dim stride
  return sts==m_scalar_type_size && (m_last_extent%vlen==0);
}

} // namespace scream

#endif // SCREAM_FIELD_ALLOC_PROP_HPP
