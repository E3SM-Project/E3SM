#include "field_alloc_prop.hpp"

namespace scream {

FieldAllocProp::FieldAllocProp (const FieldIdentifier& fid)
 : m_fid    (fid)
 , m_committed (false)
{
  // This is to have invalid initial sizes for the allocation
  m_scalar_type_size      = 0;
  m_alloc_size            = 0;
  m_alloc_value_type_size = 0;
  m_fast_index_alloc_size = 0;
}

void FieldAllocProp::commit ()
{
  if (m_committed) {
    // TODO: should we issue a warning? Error? For now, I simply do nothing
    return;
  }

  // Sanity checks: we must have requested at least one value type, and the identifier needs all dimensions set by now.
  error::runtime_check(m_alloc_value_type_size==0, "Error! No value types requested for the allocation.\n");
  error::runtime_check(m_fid.are_dimensions_set(), "Error! You need all field dimensions set before locking allocation properties.\n");

  int default_fast_index_alloc_size = m_fid.dims().back() * m_scalar_type_size;
  if (default_fast_index_alloc_size % m_alloc_value_type_size == 0) {
    // We're good! The stored value type size divides the default alloc size.
    // This is the case if, for instance, the value type is a Pack<T,N>, and N
    // divides the last dimension of the field.
    m_fast_index_alloc_size = default_fast_index_alloc_size;
  } else {
    // This is the case where there's a request for a value type whose size
    // does not divide the default alloc size. We need to pad the field.
    int num_blocks = (default_fast_index_alloc_size + m_alloc_value_type_size - 1) / m_alloc_value_type_size;
    m_fast_index_alloc_size = num_blocks * m_alloc_value_type_size;
  }

  m_alloc_size = 1;
  for (int i=0 ; i<m_fid.rank()-1; ++i) {
    m_alloc_size *= m_fid.dim(i);
  }
  m_alloc_size *= m_fast_index_alloc_size;

  m_committed = true;
}

} // namespace scream
