#include "share/field/field_identifier.hpp"

#include <ekat_string_utils.hpp>

namespace scream
{

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const Units& units,
                 const std::string& grid_name)
 : FieldIdentifier(name,layout,units,grid_name,DataType::RealType)
{
  // Nothing to do here
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const Units& units,
                 const std::string& grid_name,
                 const DataType data_type)
 : m_name      (name)
 , m_layout    (layout)
 , m_units     (units)
 , m_grid_name (grid_name)
 , m_data_type (data_type)
{
  update_identifier();
}

FieldIdentifier
FieldIdentifier::
clone (const std::string& name) const
{
  auto fid = *this;
  fid.m_name = name;
  fid.update_identifier();
  return fid;
}

FieldIdentifier& FieldIdentifier::
reset_layout (const FieldLayout& layout)
{
  m_layout = layout;
  update_identifier();
  return *this;
}

FieldIdentifier& FieldIdentifier::
reset_units  (const Units& units)
{
  m_units = units;
  update_identifier();
  return *this;
}

FieldIdentifier& FieldIdentifier::
reset_grid   (const std::string& grid)
{
  m_grid_name = grid;
  update_identifier();
  return *this;
}

FieldIdentifier& FieldIdentifier::
reset_dtype  (const DataType dtype)
{
  m_data_type = dtype;
  update_identifier();
  return *this;
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier = m_name + "[" + m_grid_name + "] <" + e2str(m_data_type);
  if (m_layout.rank()>0) {
    m_identifier += ":" + ekat::join(m_layout.names(),",");
  }
  m_identifier += ">";
  m_identifier += "(" + ekat::join(m_layout.dims(),",") + ")";

  m_identifier += " [" + m_units.to_string() + "]";
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
