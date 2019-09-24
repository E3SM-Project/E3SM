#include "field_tracking.hpp"

#include <algorithm>  // For std::find

namespace scream {

void FieldTracking::add_provider (const std::weak_ptr<AtmosphereProcess>& provider) {
  m_providers.insert(provider);
}

void FieldTracking::add_customer (const std::weak_ptr<AtmosphereProcess>& customer) {
  m_customers.insert(customer);
}

void FieldTracking::add_to_group (const std::string& group_name) {
  m_groups.insert(group_name);
}

void FieldTracking::update_time_stamp (const util::TimeStamp& ts) {
  // We check that the given time stamp is not in the past.
  // This is to prevent users from tampering with time stamps (e.g., rewinding time).
  scream_require_msg(!m_time_stamp.is_valid() || !(ts<m_time_stamp), "Error! Input time stamp is in the past.\n");

  m_time_stamp = ts;
}

} // namespace scream
