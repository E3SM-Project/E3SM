#include "share/field/field_tracking.hpp"

namespace scream {

void FieldTracking::add_provider (const std::weak_ptr<AtmosphereProcess>& provider) {
  m_providers.insert(provider);
}

void FieldTracking::add_customer (const std::weak_ptr<AtmosphereProcess>& customer) {
  m_customers.insert(customer);
}

void FieldTracking::
add_to_group (const std::shared_ptr<const FieldGroupInfo>& group) {
  m_groups.insert(group);
}

void FieldTracking::update_time_stamp (const TimeStamp& ts) {
  // We check that the given time stamp is not in the past.
  // This is to prevent users from tampering with time stamps (e.g., rewinding time).
  EKAT_REQUIRE_MSG(!m_time_stamp.is_valid() || !(ts<m_time_stamp),
      "Error! Input time stamp is in the past.\n");

  m_time_stamp = ts;

  // If you update a field, all its subviews will automatically be updated
  for (auto it : this->get_children()) {
    auto c = it.lock();
    EKAT_REQUIRE_MSG(c, "Error! A weak pointer of a child field expired.\n");
    c->update_time_stamp(ts);
  }
}

} // namespace scream
