#include "field_header.hpp"

#include "error_defs.hpp"

namespace scream
{

FieldHeader::FieldHeader (const std::string& name,
                          const std::vector<FieldTag>& tags)
 : m_name (name)
 , m_rank (tags.size())
 , m_dims (0)
 , m_tags (tags)
{
  // Nothing to be done here
}

FieldHeader::FieldHeader (const std::string& name,
                          const std::vector<int>& dims,
                          const std::vector<FieldTag>& tags)
 : FieldHeader(name,tags)
{
  set_dimensions(dims);
}

void FieldHeader::set_dimensions (const std::vector<int>& dims) {
  // We don't allow resetting dimensions after they've been set once.
  error::runtime_check(m_dims.size()==0,
                       "Error! You cannot reset field dimensions once set.\n", -1);

  // Check, then set dims
  error::runtime_check(static_cast<int>(dims.size())==m_rank,
                       "Error! Input dimensions vector not properly sized.", -1);
  for (int idim=0; idim<m_rank; ++idim) {
    error::runtime_check(dims[idim]>0, "Error! Dimensions must be positive.", -1);
  }
  m_dims = dims;
}

bool operator== (const FieldHeader& fh1, const FieldHeader& fh2) {
  // First, compare names
  if (fh1.m_name != fh2.m_name) { return false; }

  // Then, compare ranks
  if (fh1.m_rank != fh2.m_rank) { return false; }

  // Check tags along each dimension
  for (int dim=0; dim<fh1.m_rank; ++dim) {
    if (fh1.m_tags[dim]!=fh2.m_tags[dim]) { return false; }
  }

  // Check extents along each dimension
  for (int dim=0; dim<fh1.m_rank; ++dim) {
    if (fh1.m_dims[dim]!=fh2.m_dims[dim]) { return false; }
  }

  // Field are exactly the same. Return true.
  return true;
}

} // namespace scream
