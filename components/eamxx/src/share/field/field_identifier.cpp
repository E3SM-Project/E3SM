#include "share/field/field_identifier.hpp"

#include <ekat_string_utils.hpp>

namespace scream
{

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const DataType data_type)
 : FieldIdentifier(name,layout,Units::invalid(),data_type)
{
  // Nothing to do here
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const Units& units,
                 const DataType data_type)
 : m_name      (name)
 , m_layout    (layout)
 , m_units     (units)
 , m_data_type (data_type)
{
  update_identifier();
}

FieldIdentifier
FieldIdentifier::
alias (const std::string& name) const
{
  auto fid = *this;
  fid.m_name = name;
  return fid;
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier  = m_name;
  m_identifier += "<" + e2str(m_data_type) + ":" + ekat::join(m_layout.names(),",") + ">";
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
