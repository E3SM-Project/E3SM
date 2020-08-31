#include "share/field/field_header.hpp"

namespace scream
{

FieldHeader::FieldHeader (const identifier_type& id)
 : m_identifier (id)
 , m_tracking (m_identifier.name())
 , m_alloc_prop (id.get_layout())
{
  // Nothing to be done here
}

void FieldHeader::
set_extra_data (const std::string& key,
                const ekat::util::any& data,
                const bool throw_if_existing)
{
  if (throw_if_existing) {
    EKAT_REQUIRE_MSG (m_extra_data.find(key)==m_extra_data.end(),
                        "Error! Key '" + key + "' already existing in "
                        "the extra data map of field '" + m_identifier.get_id_string() + "'.\n");
    m_extra_data[key] = data;
  } else {
    m_extra_data[key] = data;
  }
}

} // namespace scream
