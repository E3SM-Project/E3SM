#include "field_tracking.hpp"

#include <algorithm>  // For std::find

namespace scream {

void FieldTracking::add_provider (const std::weak_ptr<AtmosphereProcess>& provider) {
  // Add the provider only if not already present.
  for (const auto& ptr : m_providers) {
    if (!ptr.owner_before(provider) && !provider.owner_before(ptr)) {
      return;
    }
  }
  m_providers.push_back(provider);
}

void FieldTracking::add_customer (const std::weak_ptr<AtmosphereProcess>& customer) {
  // Add the customer only if not already present.
  for (const auto& ptr : m_customers) {
    if (!ptr.owner_before(customer) && !customer.owner_before(ptr)) {
      return;
    }
  }
  m_customers.push_back(customer);
}

void FieldTracking::add_to_group (const std::string& group) {
  if (std::find(m_groups.begin(),m_groups.end(), group)==m_groups.end()) {
    m_groups.push_back(group);
  }
}

} // namespace scream
