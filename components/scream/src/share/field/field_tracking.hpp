#ifndef SCREAM_FIELD_TRACKING_HPP
#define SCREAM_FIELD_TRACKING_HPP

#include "share/util/time_stamp.hpp"

#include <memory>   // For std::weak_ptr
#include <vector>
#include <string>

namespace scream {

// Forward declarations
class AtmosphereProcess;

class FieldTracking {
public:

  // ----- Getters ----- //

  // The time stamp of the field. This can be used to check when it was last update.
  // Please, notice this is not the OS time stamp (see TimeStamp.hpp for details).
  const util::TimeStamp& get_time_stamp    () const { return m_time_stamp; }

  // List of providers/customers for this field
  const std::vector<std::weak_ptr<AtmosphereProcess>>& get_providers () const { return m_providers; }
  const std::vector<std::weak_ptr<AtmosphereProcess>>& get_customers () const { return m_customers; }

  // ----- Setters ----- //

  // Add to the list of providers/customers
  void add_provider (const std::weak_ptr<AtmosphereProcess>& provider);
  void add_customer (const std::weak_ptr<AtmosphereProcess>& customer);

  // Add the field to a given group
  void add_to_group (const std::string& group_name);

  // List of field groups that this field belongs to
  const std::vector<std::string>& get_groups_list () const { return m_groups; }

protected:
  

  // Tracking the updates of the field
  util::TimeStamp   m_time_stamp;

  // These are to be used to track the order in which providers update the field at each time step.
  // One can use this information to track when a field gets updated during a timestep. It can be
  // particularly useful in the case of parallel schedules.
  std::vector<std::string> m_last_ts_providers;
  std::vector<std::string> m_curr_ts_providers;

  // List of provider/customer processes. A provider is an atm process that computes/updates the field.
  // A customer is an atm process that uses the field just as an input.
  // NOTE: do NOT use shared_ptr, since you will likely create circular references.
  std::vector<std::weak_ptr<AtmosphereProcess>> m_providers;
  std::vector<std::weak_ptr<AtmosphereProcess>> m_customers;

  // Groups are used to bundle together fields, so that a process can request all of them
  // without knowing/listing all their names. For instance, the dynamics process needs to
  // get all tracers, which need to be advected. However, dyamics has no idea of what are
  // the tracers names, and neither should it care. Groups can come to rescue here, allowing
  // dynamics to request all fields that have been marked as 'tracers'.
  std::vector<std::string>    m_groups;
};

} // namespace scream

#endif // SCREAM_FIELD_TRACKING_HPP
