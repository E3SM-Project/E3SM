#include "share/field/field_tracking.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include <algorithm>  // For std::find

namespace scream {

FieldTracking::FieldTracking (const std::string& name)
 : m_name (name)
{
  // Nothing to do here
}

FieldTracking::FieldTracking(const std::string& name,
                             const std::shared_ptr<FieldTracking>& parent)
 : FieldTracking(name)
{
  m_parent = parent;

  // Lock the parent, and see if there are any tracking properties that should be propagated

  auto p = m_parent.lock();
  EKAT_REQUIRE_MSG (p!=nullptr, "Error! Input parent pointer not valid.\n");

  // If the parent already has a valid time stamp, set it here
  if (p->get_time_stamp().is_valid()) {
    update_time_stamp(p->get_time_stamp());
  }
}

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
  for (auto c : m_children) {
    c.lock()->update_time_stamp(ts);
  }
}

void FieldTracking::register_as_children_in_parent () {
  if (m_parent.lock()==nullptr) {
    return;
  }

  // Scan the children of my parent. If I'm not already there, add myself.
  auto me = shared_from_this();
  auto siblings = m_parent.lock()->m_children;
  for (auto p : siblings) {
    if (p.lock()==me) {
      return;
    }
  }

  siblings.push_back(me);
}

} // namespace scream
