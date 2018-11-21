#include "field_alloc_prop.hpp"

namespace scream {

FieldAllocProp::FieldAllocProp (const FieldIdentifier& fid)
 : m_fid    (fid)
 , m_committed (false)
{
  // This is to have invalid initial sizes for the allocation
  m_scalar_type_name = "";
  m_scalar_type_size = 0;
  m_alloc_size       = 0;
  m_value_type_size  = 0;
}

void FieldAllocProp::commit ()
{
  if (m_committed) {
    // TODO: should we issue a warning? Error? For now, I simply do nothing
    return;
  }

  // Sanity checks: we must have requested at least one value type, and the identifier needs all dimensions set by now.
  error::runtime_check(m_value_type_size>0, "Error! No value types requested for the allocation.\n");
  error::runtime_check(m_fid.are_dimensions_set(), "Error! You need all field dimensions set before locking allocation properties.\n");

  int last_dim = m_fid.dims().back();
  int default_last_dim_alloc_size = last_dim*m_scalar_type_size;
  if (default_last_dim_alloc_size % m_value_type_size == 0) {
    // We're good! The stored value type size divides the default alloc size.
    // This is the case if, for instance, the value type is a Pack<T,N>, and N
    // divides the last dimension of the field.
  } else {
    // This is the case where there's a request for a value type whose size
    // does not divide the default alloc size. We need to pad the field.
    last_dim = (default_last_dim_alloc_size + m_value_type_size - 1) / m_value_type_size;
    // last_dim_alloc_size = last_dim * m_value_type_size;
  }

  // Determining the total allocation
  m_alloc_size = m_value_type_size;;
  for (int i=0 ; i<m_fid.rank()-1; ++i) {
    m_alloc_size *= m_fid.dim(i);
  }
  m_alloc_size *= last_dim;

  m_committed = true;
}

} // namespace scream
