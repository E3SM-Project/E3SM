#include "field_alloc_prop.hpp"

namespace scream {

FieldAllocProp::FieldAllocProp (const FieldIdentifier& fid)
 : m_fid       (fid)
 , m_committed (false)
{
  // This is to have invalid initial sizes for the allocation
  m_scalar_type_name = "";
  m_scalar_type_size = 0;
  m_alloc_size       = 0;
}

void FieldAllocProp::commit ()
{
  if (m_committed) {
    // TODO: should we issue a warning? Error? For now, I simply do nothing
    return;
  }

  // Sanity checks: we must have requested at least one value type, and the identifier needs all dimensions set by now.
  error::runtime_check(m_value_type_sizes.size()>0, "Error! No value types requested for the allocation.\n");
  error::runtime_check(m_fid.are_dimensions_set(), "Error! You need all field dimensions set before committing the allocation properties.\n");

  // Loop on all value type sizes.
  for (auto vts : m_value_type_sizes) {
    // The number of scalar_type in a value_type
    const int num_scalars_in_value_type = vts / m_scalar_type_size;

    // The number of value_type's needed in the fast-striding dimension
    const int num_last_dim_value_types = (m_fid.dims().back() + num_scalars_in_value_type - 1) / num_scalars_in_value_type;

    // The size of such allocation
    const int value_type_alloc_size = (m_fid.size() / m_fid.dims().back()) // All except the last dimension
                                    * num_last_dim_value_types*vts;

    m_alloc_size = std::max(m_alloc_size, value_type_alloc_size);
  }

  m_committed = true;
}

} // namespace scream
