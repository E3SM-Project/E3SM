#include "field_alloc_prop.hpp"

namespace scream {

FieldAllocProp::FieldAllocProp ()
 : m_value_type_sizes (0)
 , m_scalar_type_size (0)
 , m_scalar_type_name ("")
 , m_alloc_size       (0)
 , m_committed        (false)
{
  // Nothing to be done here
}

void FieldAllocProp::commit (const FieldLayout& layout)
{
  if (m_committed) {
    // TODO: should we issue a warning? Error? For now, I simply do nothing
    return;
  }

  // Sanity checks: we must have requested at least one value type, and the identifier needs all dimensions set by now.
  EKAT_REQUIRE_MSG(m_value_type_sizes.size()>0, "Error! No value types requested for the allocation.\n");
  EKAT_REQUIRE_MSG(layout.are_dimensions_set(), "Error! You need all field dimensions set before committing the allocation properties.\n");

  // Loop on all value type sizes.
  m_last_dim_alloc_size = 0;
  for (auto vts : m_value_type_sizes) {
    // The number of scalar_type in a value_type
    const int num_scalars_in_value_type = vts / m_scalar_type_size;

    // The number of value_type's needed in the fast-striding dimension
    const int num_last_dim_value_types = (layout.dims().back() + num_scalars_in_value_type - 1) / num_scalars_in_value_type;

    // The size of such allocation (along the last dim only)
    m_last_dim_alloc_size = std::max(m_last_dim_alloc_size, vts * num_last_dim_value_types);
  }

  m_alloc_size = (layout.size() / layout.dims().back()) // All except the last dimension
                 * m_last_dim_alloc_size;

  m_committed = true;
}

} // namespace scream
