#include "field_header.hpp"

#include "error_defs.hpp"

namespace scream
{

FieldHeader::FieldHeader (const std::string& name,
                          const std::vector<FieldTag>& tags)
 : m_name (name)
 , m_tags (tags)
 , m_dims (0)
 , m_rank (tags.size())
{
  update_identifier ();
}

FieldHeader::FieldHeader (const std::string& name,
                          const std::vector<FieldTag>& tags,
                          const std::vector<int>& dims)
 : FieldHeader(name,tags)
{
  set_dimensions(dims);
  update_identifier ();
}

void FieldHeader::set_dimensions (const std::vector<int>& dims) {
  // We don't allow resetting dimensions after they've been set once.
  error::runtime_check(!dimensions_set(),
                       "Error! You cannot reset field dimensions once set.\n", -1);

  // Check, then set dims
  error::runtime_check(dims.size()==static_cast<size_t>(m_rank),
                       "Error! Input dimensions vector not properly sized.", -1);
  for (int idim=0; idim<m_rank; ++idim) {
    error::runtime_check(dims[idim]>0, "Error! Dimensions must be positive.", -1);
  }
  m_dims = dims;
}

void FieldHeader::add_provider (std::shared_ptr<AtmosphereProcess> provider) {
  m_providers.insert(provider);
}

void FieldHeader::add_customer (std::shared_ptr<AtmosphereProcess> customer) {
  m_customers.insert(customer);
}

void FieldHeader::update_identifier () {
  // Update the identifier
  m_identifier = m_name + "<" + tag2string(m_tags[0]);
  for (int dim=1; dim<m_rank; ++dim) {
    m_identifier += "," + tag2string(m_tags[dim]);
  }
  m_identifier += "> (";
  if (m_dims.size()>0) {
    m_identifier += m_dims[0];
    for (int dim=1; dim<m_rank; ++dim) {
      m_identifier += "," + std::to_string(m_dims[dim]);
    }
  }
  m_identifier += ")";
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
