#include "share/field/field_tracking.hpp"
#include "share/field/field_initializer.hpp"
#include "share/field/field_value_initializer.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include <algorithm>  // For std::find

namespace scream {

FieldTracking::FieldTracking (const std::string& name)
 : m_name (name)
 , m_init_type (InitType::None)
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

  // If an initializer was already set in the parent, make sure
  // this field doesn't get assigned another initializer
  if (p->get_init_type()!=InitType::None) {
    set_init_type(InitType::Inherited);
  }

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

void FieldTracking::set_initializer (const std::weak_ptr<FieldInitializer>& initializer) {
  // Try to set the init type. If init type was already set, or not needed, an error will be raised.
  set_init_type (InitType::Initializer);

  // Make sure nobody else already claimed this role
  EKAT_REQUIRE_MSG (!static_cast<bool>(m_initializer.lock()),
      "Error! There is already an initializer for field '" + m_name + "'.\n" +
      "       Current initializer: " + m_initializer.lock()->name() + "\n");

  // Make sure the initializer is nonnull
  EKAT_REQUIRE_MSG (static_cast<bool>(initializer.lock()),
      "Error! Input initializer process is invalid.\n");

  m_initializer = initializer;

  // If you init a field, all its subviews will automatically be inited
  for (auto c : m_children) {
    // Can't call set_initializer on c, since it would set init type to Initializer
    // All we need is to mark that the field is somehow inited.
    c.lock()->set_init_type(InitType::Inherited);
  }
}

void FieldTracking::set_value_initializer (const Real value) {

  // Try to set the init type. If init type was already set, or not needed, an error will be raised.
  set_init_type (InitType::Value);

  // Make sure nobody else already claimed this role
  EKAT_REQUIRE_MSG (!static_cast<bool>(m_initializer.lock()),
      "Error! There is already an initializer for field '" + m_name + "'.\n" +
      "       Current initializer: " + m_initializer.lock()->name() + "\n");

  m_initializer = create_field_value_initializer(value);

  // If you init a field, all its subviews will automatically be inited
  for (auto c : m_children) {
    // Can't call set_initializer on c, since it would set init type to Initializer
    // All we need is to mark that the field is somehow inited.
    c.lock()->set_init_type(InitType::Inherited);
  }
}

void FieldTracking::set_init_type (const InitType init_type) {

  // EKAT_REQUIRE_MSG(m_init_type!=InitType::NotNeeded,
  //   "Error! This field was not expected to need an initialization.\n");

  // The following changes are the only allowed ones
  //  - None -> X = X
  //  - X -> X = X
  //  - X -> None = X
  // Any other attempt to modify the init type will result in an error.
  // In particular, if m_init_type is NotNeeded, calling this method is an error.
  if (m_init_type==InitType::None) {
    m_init_type = init_type;
  } else if (init_type!=InitType::None) {
    EKAT_REQUIRE_MSG (init_type==m_init_type,
        "Error! Invalid attempt made to modify initialization type for field '" + m_name + "'.\n" +
        "       Old init type: " + e2str(m_init_type) + "\n" +
        "       New init type: " + e2str(init_type) + "\n");
  }
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
