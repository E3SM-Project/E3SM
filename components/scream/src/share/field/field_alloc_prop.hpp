#ifndef SCREAM_FIELD_ALLOC_PROP_HPP
#define SCREAM_FIELD_ALLOC_PROP_HPP

#include "share/field/field_identifier.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/ekat_assert.hpp"

#include <vector>

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
 *      In the field repo, this would be stored as a 1d field: View<Real*>.
 *      Let's say you have one class that uses the field as View<Real**>,
 *      and another class that used the field in a packed way, taht is,
 *      as View<Pack<Real,4>**>. The second view will have dimensions (3,4),
        in terms of its value type (i.e., Pack<Real,4>).
 *      This means the field needs an allocation bigger than 30 Real's, namely
 *      3*4*4=36 Real's. This class does the book-keeping for the allocation size,
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
 *      users need to call the get_reshaped_view<T>() to use it with
 *      a different data type. E.g., to use the field as a 2d field,
 *      with Pack<Real,8> as scalar type, one needs to call
 *        auto v = f.get_reshaped_view<Pack<Real,8>**>()
 *      The alloc props are queried to establish 1) whether the allocation
 *      is compatible with a pack-8 scalar type, and 2) how the dimensions
 *      of the resulting view should be.
 *    - It is possible (with some limitations) to have a Field being a
 *      "subview" of a bigger field. E.g., tracers are allocated as a big
 *      array Q, but one may just be interested in a single one (e.g., qv),
 *      yet still use it as a Field. This requires some care during the
 *      get_reshaped_view method, so the field alloc props can be used
 *      to detect whether we are in this scenario.
 *
 *  Note: at every request for a new value_type, this class checks the
 *        underlying scalar_type. We ASSUME that the dimensions in the
 *        field identifier refer to that scalar_type.
 */

class FieldAllocProp {
public:
  using layout_type      = FieldLayout;
  using layout_ptr_type  = std::shared_ptr<const layout_type>;

  FieldAllocProp ();

  FieldAllocProp& operator= (const FieldAllocProp&);

  // Return allocation props of an unmanaged subivew of this field
  // at entry k along dimension idim.
  FieldAllocProp subview (const int idim, const int k) const;

  // Request allocation able to accommodate the given ValueType
  template<typename ValueType>
  void request_allocation ();

  // Request allocation able to accommodate a pack of ScalarType of the given pack size
  template<typename ScalarType>
  void request_allocation (int pack_size);

  // Request allocation able to accommodate pack_size scalars,
  // where scalars have size scalar_size
  void request_allocation (int scalar_size, int pack_size);

  // Request allocation able to accommodate all the alloc props in src.
  // Note: src does not need to be committed yet.
  void request_allocation (const FieldAllocProp& src);

  // Locks the allocation properties, preventing furter value types requests
  void commit (const layout_ptr_type& layout);

  // ---- Getters ---- //

  // Whether commit has been called
  bool is_committed () const { return m_committed; }

  // Get the overall allocation size (in Bytes)
  int  get_alloc_size () const;

  // Wether this allocation is contiguous
  bool contiguous () const { return m_contiguous; }

  // Size of the last extent in the alloction (i.e., number of scalars in it)
  int  get_last_extent () const;
  int  get_padding () const;

  // Return the slice information of this subview allocation.
  const std::pair<int,int>& get_subview_idx () const { return m_subview_idx; }

  // Whether the allocation properties are compatible with ValueType
  template<typename ValueType>
  bool is_compatible () const;

protected:

  // The FieldLayout associated to this allocation
  layout_ptr_type     m_layout;

  // The list of requested value types for this allocation
  std::vector<int>    m_value_type_sizes;

  // The size of the scalar type
  int   m_scalar_type_size;

  // The extent along the last dimension for this allocation
  int   m_last_extent;

  // The total size of this allocation.
  int   m_alloc_size;

  // If this allocation is a subview of another, this pair stores the
  // dimension idim and entry k along dimension idim of the slice.
  std::pair<int,int>  m_subview_idx;

  // Whether the allocation is contiguous
  bool  m_contiguous;

  // Whether commit was called (i.e., no more value type requests allowed)
  bool  m_committed;
};

// ================================= IMPLEMENTATION ================================== //

template<typename ValueType>
void FieldAllocProp::request_allocation () {
  using ekat::ScalarTraits;
  using scalar_type = typename ScalarTraits<ValueType>::scalar_type;
  const int n = sizeof(ValueType) / sizeof(scalar_type);
  request_allocation<scalar_type>(n);
}

template<typename ScalarType>
void FieldAllocProp::request_allocation (const int pack_size) {
  request_allocation (sizeof(ScalarType),pack_size);
}

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
