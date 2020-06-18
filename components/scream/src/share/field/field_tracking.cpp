#include "share/field/field_tracking.hpp"
#include "share/field/field_initializer.hpp"
#include "share/field/field_value_initializer.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include <algorithm>  // For std::find

namespace scream {

FieldTracking::FieldTracking (const std::string& name)
 : m_name (name)
 // , m_init_type (InitType::NotNeeded)
 , m_init_type (InitType::None)
{
  // Nothing else to do
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
  scream_require_msg (!static_cast<bool>(m_initializer.lock()),
                      "Error! There is already an initializer for field '" + m_name + "'.\n" +
                      "       Current initializer: " + m_initializer.lock()->name() + "\n");

  // Make sure the initializer is nonnull
  scream_require_msg (static_cast<bool>(initializer.lock()),
                      "Error! Input initializer process is invalid.\n");

  m_initializer = initializer;
}

void FieldTracking::set_value_initializer (const Real value) {

  // Try to set the init type. If init type was already set, or not needed, an error will be raised.
  set_init_type (InitType::Value);

  // Make sure nobody else already claimed this role
  scream_require_msg (!static_cast<bool>(m_initializer.lock()),
                      "Error! There is already an initializer for field '" + m_name + "'.\n" +
                      "       Current initializer: " + m_initializer.lock()->name() + "\n");

  m_initializer = create_field_value_initializer(value);
}

void FieldTracking::set_init_type (const InitType init_type) {

  // scream_require_msg(m_init_type!=InitType::NotNeeded,
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
    scream_require_msg (init_type==m_init_type,
                        "Error! Invalid attempt made to modify initialization type for field '" + m_name + "'.\n" +
                        "       Old init type: " + e2str(m_init_type) + "\n" +
                        "       New init type: " + e2str(init_type) + "\n");
  }
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
