#include "field_identifier.hpp"

namespace scream
{

FieldIdentifier::FieldIdentifier (const std::string& name,
                                  const layout_type& layout)
 : m_name   (name)
 , m_layout (layout)
{
  update_identifier ();
}

FieldIdentifier::FieldIdentifier (const std::string& name,
                                  const std::vector<FieldTag>& tags)
 : m_name   (name)
 , m_layout (tags)
{
  // Compute identifier string
  update_identifier ();
}

FieldIdentifier::FieldIdentifier (const std::string& name,
                                  const std::initializer_list<FieldTag>& tags)
 : m_name   (name)
 , m_layout (tags)
{
  update_identifier ();
}

void FieldIdentifier::set_dimension (const int idim, const int dimension) {
  m_layout.set_dimension(idim,dimension);
  update_identifier ();
}

void FieldIdentifier::set_dimensions (const std::vector<int>& dims) {
  m_layout.set_dimensions(dims);
  update_identifier ();
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier = m_name + "<" + tag2string(m_layout.tags()[0]);
  for (int dim=1; dim<m_layout.rank(); ++dim) {
    m_identifier += "," + tag2string(m_layout.tags()[dim]);
  }
  m_identifier += ">(" + std::to_string(m_layout.dims()[0]);
  for (int dim=1; dim<m_layout.rank(); ++dim) {
    m_identifier += "," + std::to_string(m_layout.dims()[dim]);
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
