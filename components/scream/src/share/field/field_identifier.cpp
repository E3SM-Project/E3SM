#include "share/field/field_identifier.hpp"

namespace scream
{

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const Units& units,
                 const std::string& grid_name)
 : m_name   (name)
 , m_layout (layout)
 , m_units  (units)
{
  // This also calls 'update_identifier'
  set_grid_name(grid_name);
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const std::vector<FieldTag>& tags,
                 const Units& units,
                 const std::string& grid_name)
 : m_name   (name)
 , m_layout (tags)
 , m_units  (units)
{
  // This also calls 'update_identifier'
  set_grid_name(grid_name);
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const std::initializer_list<FieldTag>& tags,
                 const Units& units,
                 const std::string& grid_name)
 : m_name   (name)
 , m_layout (tags)
 , m_units  (units)
{
  // This also calls 'update_identifier'
  set_grid_name(grid_name);
}

void FieldIdentifier::set_dimension (const int idim, const int dimension) {
  m_layout.set_dimension(idim,dimension);
  update_identifier ();
}

void FieldIdentifier::set_dimensions (const std::vector<int>& dims) {
  m_layout.set_dimensions(dims);
  update_identifier ();
}

void FieldIdentifier::set_grid_name (const std::string& grid_name) {
  // Only allow overwriting if the stored grid name is empty
  EKAT_REQUIRE_MSG (m_grid_name=="", "Error! Cannot overwrite a non-empty grid name.\n");

  m_grid_name = grid_name;

  // Update the identifier string
  update_identifier ();
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier = m_name + "[" + m_grid_name + "]<" + tag2string(m_layout.tags()[0]);
  for (int dim=1; dim<m_layout.rank(); ++dim) {
    m_identifier += "," + tag2string(m_layout.tags()[dim]);
  }
  m_identifier += ">(" + std::to_string(m_layout.dims()[0]);
  for (int dim=1; dim<m_layout.rank(); ++dim) {
    m_identifier += "," + std::to_string(m_layout.dims()[dim]);
  }
  m_identifier += ") [" + m_units.get_string() + "]";
}

// Free functions for identifiers comparison
bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_id_string()==fid2.get_id_string());
}

bool operator< (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_id_string()<fid2.get_id_string());
}

} // namespace scream
