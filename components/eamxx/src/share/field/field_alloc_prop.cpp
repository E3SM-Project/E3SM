#include "field_alloc_prop.hpp"

namespace scream {

FieldAllocProp::FieldAllocProp (const int scalar_size)
 : m_value_type_sizes (1,scalar_size)
 , m_scalar_type_size (scalar_size)
 , m_pack_size_max    (1)
 , m_alloc_size       (0)
 , m_committed        (false)
{
  // Nothing to do here
}

FieldAllocProp& FieldAllocProp::operator= (const FieldAllocProp& src)
{
  if  (&src!=this) {
    EKAT_REQUIRE_MSG (not m_committed,
        "Error! Cannot assign FieldAllocProp once the dst obj is committed.\n");

    m_layout = src.m_layout;
    m_value_type_sizes = src.m_value_type_sizes;
    m_scalar_type_size = src.m_scalar_type_size;
    m_pack_size_max = src.m_pack_size_max;
    m_last_extent = src.m_last_extent;
    m_alloc_size = src.m_alloc_size;
    m_subview_info = src.m_subview_info;
    m_contiguous = src.m_contiguous;
    m_committed = src.m_committed;
  }
  return *this;
}

FieldAllocProp FieldAllocProp::
subview (const int idim, const int k, const bool dynamic) const {
  EKAT_REQUIRE_MSG(is_committed(),
      "Error! Subview requires alloc properties to be committed.\n");
  EKAT_REQUIRE_MSG (idim==0 || idim==1,
      "Error! Subviewing is only allowed along first or second dimension.\n");
  EKAT_REQUIRE_MSG (idim<m_layout->rank(),
      "Error! Dimension index out of bounds.\n");
  EKAT_REQUIRE_MSG (k>=0 && k<m_layout->dim(idim),
      "Error! Index along the dimension is out of bounds.\n");

  // Set new layout basic stuff
  FieldAllocProp props(m_scalar_type_size);
  props.m_committed = false;
  props.m_scalar_type_size = m_scalar_type_size;
  props.m_pack_size_max = m_pack_size_max;
  props.m_alloc_size = m_alloc_size / m_layout->dim(idim);
  props.m_subview_info = SubviewInfo(idim,k,m_layout->dim(idim),dynamic);

  // The output props should still store a FieldLayout, in case
  // they are further subviewed. We have all we need here to build
  // a layout from scratch, but it would duplicate what is inside
  // the FieldIdentifier that will be built for the corresponding
  // field. Therefore, we'll require the user to still call 'commit',
  // passing the layout from the field id.
  // HOWEVER, this may cause bugs, cause the user might make a mistake,
  // and pass the wrong layout. Therefore, we build a "temporary" one,
  // which will be replaced during the call to 'commit'.
  std::vector<FieldTag> tags = m_layout->tags();
  std::vector<int> dims = m_layout->dims();
  tags.erase(tags.begin()+idim);
  dims.erase(dims.begin()+idim);
  props.m_layout = std::make_shared<layout_type>(tags,dims);

  // Output is contioguous if either
  //  - this->m_contiguous=true AND idim==0
  //  - m_layout->dim(i)==1 for all i<idim
  //  - props.m_layout->rank()==0
  props.m_contiguous = m_contiguous || props.m_layout->rank()==0;
  for (int i=0; i<idim; ++i) {
    if (m_layout->dim(i)>0) {
      props.m_contiguous = false;
      break;
    }
  }

  // Figure out strides
  const int rm1 = m_layout->rank()-1;
  if (idim==rm1) {
    // We're slicing the possibly padded dim, so everything else is as in the layout
    props.m_last_extent = m_layout->dim(idim);
  } else {
    // We are keeping the last dim, so same last extent
    props.m_last_extent = m_last_extent;
  }
  return props;
}

void FieldAllocProp::request_allocation (const int pack_size) {
  using ekat::ScalarTraits;

  EKAT_REQUIRE_MSG(!m_committed,
      "Error! Cannot change allocation properties after they have been commited.\n");

  const int vts = m_scalar_type_size*pack_size;

  // Store the size of the value type.
  m_value_type_sizes.push_back(vts);
}

void FieldAllocProp::request_allocation (const FieldAllocProp& src)
{
  const auto sts = src.m_scalar_type_size;
  for (auto vts : src.m_value_type_sizes) {
    // For each value type size, the pack size is simply vts/sts.
    request_allocation(vts/sts);
  }
}

int FieldAllocProp::get_padding () const {
  EKAT_REQUIRE_MSG(is_committed(),
      "Error! You cannot query the allocation padding until after calling commit().");
  int padding = m_last_extent - m_layout->dims().back();
  return padding;
}

void FieldAllocProp::reset_subview_idx (const int idx) {
  EKAT_REQUIRE_MSG (is_committed(),
      "Error! Cannot reset subview idx on a non-committed allocation.\n");
  EKAT_REQUIRE_MSG (is_subfield(),
      "Error! Cannot reset subview idx if this is not a subfield.\n");
  EKAT_REQUIRE_MSG (is_dynamic_subfield(),
      "Error! Cannot reset subview idx for non-dynamic subfields.\n");
  EKAT_REQUIRE_MSG (idx>=0 && idx<m_subview_info.dim_extent,
      "Error! Subview slice idx out of bounds.\n");

  // Note: all the other allocation properties are unchanged.
  m_subview_info.slice_idx = idx;
}

void FieldAllocProp::commit (const layout_ptr_type& layout)
{
  if (is_committed()) {
    // TODO: should we issue a warning? Error? For now, I simply do nothing
    return;
  }

  EKAT_REQUIRE_MSG (layout,
      "Error! Invalid input layout pointer.\n");

  if (m_alloc_size>0) {
    // This obj was created as a subview of another alloc props obj.
    // Check that input layout matches the stored one, then replace the ptr.
    EKAT_REQUIRE_MSG (*layout==*m_layout,
        "Error! The input field layout does not match the stored one.\n");

    m_layout = layout;
    m_committed = true;
    return;
  }

  // Sanity checks: we must have requested at least one value type, and the identifier needs all dimensions set by now.
  EKAT_REQUIRE_MSG(m_value_type_sizes.size()>0,
      "Error! No value types requested for the allocation.\n");
  EKAT_REQUIRE_MSG(layout->are_dimensions_set(),
      "Error! You need all field dimensions set before committing the allocation properties.\n");

  // Store pointer to layout for future use (in case subview is called)
  m_layout = layout;

  // Loop on all value type sizes.
  m_last_extent = 0;
  int last_phys_extent = m_layout->dims().back();
  for (auto vts : m_value_type_sizes) {
    // The number of scalar_type in a value_type
    const int vt_len = vts / m_scalar_type_size;

    // Update the max pack size
    m_pack_size_max = std::max(m_pack_size_max,vt_len);

    // The number of value_type's needed in the fast-striding dimension
    const int num_vt = (last_phys_extent + vt_len - 1) / vt_len;

    // The total number of scalars with num_vt value_type entries
    const int num_st = num_vt*vt_len;

    // The size of such allocation (along the last dim only)
    m_last_extent = std::max(m_last_extent, num_st);
  }

  // If we have a partitioned grid, with some ranks not owning any grid point,
  // we may end up with a layout of size 0
  if (m_layout->size()==0) {
    m_alloc_size = 0;
  } else {
    m_alloc_size = (m_layout->size() / last_phys_extent) // All except the last dimension
                   * m_last_extent * m_scalar_type_size;
  }

  m_contiguous = true;

  m_committed = true;
}

long long FieldAllocProp::get_alloc_size () const {
  EKAT_REQUIRE_MSG(is_committed(),
      "Error! You cannot query the allocation properties until they have been committed.");
  return m_alloc_size;
}

int FieldAllocProp::get_largest_pack_size () const {
  EKAT_REQUIRE_MSG(is_committed(),
      "Error! You cannot query the allocation properties until they have been committed.");
  return m_pack_size_max;
}

int FieldAllocProp::get_last_extent () const {
  EKAT_REQUIRE_MSG(is_committed(),
      "Error! You cannot query the allocation strides until after calling commit().");
  return m_last_extent;
}

} // namespace scream
