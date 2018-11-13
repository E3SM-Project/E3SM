#include "field_header.hpp"

#include <share/error_defs.hpp>
#include <algorithm>  // For std::find

namespace scream
{

FieldHeader::FieldHeader (const identifier_type& id)
 : m_identifier(id)
{
  // Nothing to be done here
}

void FieldHeader::add_provider (std::shared_ptr<AtmosphereProcess> provider) {
  if (std::find(m_providers.begin(),m_providers.end(), provider)==m_providers.end()) {
    m_providers.push_back(provider);
  }
}

void FieldHeader::add_customer (std::shared_ptr<AtmosphereProcess> customer) {
  if (std::find(m_customers.begin(),m_customers.end(), customer)==m_customers.end()) {
    m_customers.push_back(customer);
  }
}

void FieldHeader::add_to_group (const std::string& group) {
  if (std::find(m_groups.begin(),m_groups.end(), group)==m_groups.end()) {
    m_groups.push_back(group);
  }
}

// bool operator== (const FieldHeader& fh1, const FieldHeader& fh2) {
//   // First, compare identifiers
//   if (fh1.m_name != fh2.m_name) { return false; }

//   // Then, compare ranks
//   if (fh1.m_rank != fh2.m_rank) { return false; }

//   // Check tags along each dimension
//   for (int dim=0; dim<fh1.m_rank; ++dim) {
//     if (fh1.m_tags[dim]!=fh2.m_tags[dim]) { return false; }
//   }

//   // Check extents along each dimension
//   for (int dim=0; dim<fh1.m_rank; ++dim) {
//     if (fh1.m_dims[dim]!=fh2.m_dims[dim]) { return false; }
//   }

//   // Field are exactly the same. Return true.
//   return true;
// }

} // namespace scream
