#include "share/field/field_header.hpp"

namespace scream
{

FieldHeader::FieldHeader (const identifier_type& id)
 : m_identifier (id)
 , m_alloc_prop (id.get_layout())
{
  // Nothing to be done here
}

} // namespace scream
