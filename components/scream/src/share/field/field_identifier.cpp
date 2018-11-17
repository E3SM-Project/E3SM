#include "field_identifier.hpp"

namespace scream
{

FieldIdentifier::FieldIdentifier (const std::string& name,
                                  const std::vector<FieldTag>& tags)
 : m_name (name)
 , m_tags (tags)
 , m_rank (tags.size())
 , m_dims (tags.size(),-1)
{
  // Compute identifier string
  update_identifier ();
}

void FieldIdentifier::set_dimension (const int idim, const int dimension) {
  error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  error::runtime_check(m_dims[idim] == -1,
                       "Error! You cannot reset field dimensions once set.\n");
  error::runtime_check(dimension>0, "Error! Dimensions must be positive.");
  m_dims[idim] = dimension;

  update_identifier ();
}

void FieldIdentifier::set_dimensions (const std::vector<int>& dims) {
  // Check, then set dims
  error::runtime_check(dims.size()==static_cast<size_t>(m_rank),
                       "Error! Input dimensions vector not properly sized.");
  for (int idim=0; idim<m_rank; ++idim) {
    set_dimension(idim,dims[idim]);
  }
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier = m_name + "<" + tag2string(m_tags[0]);
  for (int dim=1; dim<m_rank; ++dim) {
    m_identifier += "," + tag2string(m_tags[dim]);
  }
  m_identifier += ">(" + std::to_string(m_dims[0]);
  for (int dim=1; dim<m_rank; ++dim) {
    m_identifier += "," + std::to_string(m_dims[dim]);
  }
  m_identifier += ")";
}

// Free functions for identifiers comparison
bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_identifier()==fid2.get_identifier());
}

bool operator< (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_identifier()<fid2.get_identifier());
}

} // namespace scream
