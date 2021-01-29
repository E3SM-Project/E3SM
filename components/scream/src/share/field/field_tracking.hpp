#ifndef SCREAM_FIELD_TRACKING_HPP
#define SCREAM_FIELD_TRACKING_HPP

#include "share/scream_types.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <memory>   // For std::weak_ptr
#include <string>
#include <list>
#include <set>

namespace scream {

// How a field has to be init-ed
// Note: internally, Value still creates a FieldInitializer object, but can be done
//       behind the scenes by the infrastructure
enum class InitType {
  // NotNeeded,    // No initialization is needed for this field
  Value,        // Should be inited to a specific value
  Initializer,  // A FieldInitializer object should take care of this
  Inherited,    // For subviews: an init type has been set for the 'parent' field
  None          // No initialization is needed/expected
};

inline std::string e2str (const InitType e) {
  std::string s;
  switch (e) {
    // case InitType::NotNeeded:   s = "NotNeeded";    break;
    case InitType::None:        s = "None";         break;
    case InitType::Value:       s = "Value";        break;
    case InitType::Initializer: s = "Initializer";  break;
    case InitType::Inherited:   s = "Inherited";    break;
    default: s = "INVALID";
  }

  return s;
}

// Forward declarations
class AtmosphereProcess;
class FieldInitializer;
class FieldHeader;

class FieldTracking : public ekat::enable_shared_from_this<FieldTracking> {
public:

  using TimeStamp         = util::TimeStamp;
  using ci_string         = ekat::CaseInsensitiveString;
  using atm_proc_ptr_type = std::weak_ptr<AtmosphereProcess>;
  using atm_proc_set_type = ekat::WeakPtrSet<AtmosphereProcess>;

  FieldTracking() = delete;
  FieldTracking(const std::string& name);
  FieldTracking(const FieldTracking&) = default;
  FieldTracking(const std::string& name,
                const std::shared_ptr<FieldTracking>& parent);

  // No assignment, to prevent tampering with tracking (e.g., rewinding time stamps)
  FieldTracking& operator=(const FieldTracking&) = delete;

  // ----- Getters ----- //

  // The time stamp of the field. This can be used to check when it was last updated.
  // Please, notice this is not the OS time stamp (see time_stamp.hpp for details).
  const TimeStamp& get_time_stamp () const { return m_time_stamp; }

  // List of providers/customers for this field
  const atm_proc_set_type& get_providers () const { return m_providers; }
  const atm_proc_set_type& get_customers () const { return m_customers; }
  const std::weak_ptr<FieldInitializer>& get_initializer () const { return m_initializer; }
  InitType get_init_type () const { return m_init_type; }

  // Get parent/children (if any)
  // std::shared_ptr<const FieldHeader> get_parent () const { return m_parent.lock(); }
  std::weak_ptr<FieldTracking> get_parent () const { return m_parent; }

  // List of field groups that this field belongs to
  const std::set<ci_string>& get_groups_names () const { return m_groups; }

  // ----- Setters ----- //

  // Add to the list of providers/customers
  void add_provider (const std::weak_ptr<AtmosphereProcess>& provider);
  void add_customer (const std::weak_ptr<AtmosphereProcess>& customer);
  void set_initializer (const std::weak_ptr<FieldInitializer>& initializer);
  void set_value_initializer (const Real value);

  // Add the field to a given group
  void add_to_group (const std::string& group_name);

  // Set the time stamp for this field. This can only be called once, due to TimeStamp implementation.
  // NOTE: if the field has 'children' (see below), their ts will be udpated too.
  //       However, if the field has a 'parent' (see below), the parent's ts will not be updated.
  void update_time_stamp (const TimeStamp& ts);

  const std::string& name () const { return m_name; }

  void register_as_children_in_parent ();

protected:

  void set_init_type (const InitType init_type);

  // We keep the field name just to make debugging messages more helpful
  std::string m_name;

  // Tracking the updates of the field
  TimeStamp         m_time_stamp;

  // These are to be used to track the order in which providers update the field at each time step.
  // One can use this information to track when a field gets updated during a timestep. It can be
  // particularly useful in the case of parallel schedules.
  std::set<std::string>   m_last_timestep_providers;
  std::set<std::string>   m_curr_timestep_providers;

  // List of provider/customer processes. A provider is an atm process that computes/updates the field.
  // A customer is an atm process that uses the field just as an input.
  // NOTE: do NOT use shared_ptr, since you would create circular references.
  atm_proc_set_type       m_providers;
  atm_proc_set_type       m_customers;

  // How this field will be initialized (if at all needed)
  InitType                m_init_type;

  // The initializer is a class that claims the responsibility of initializing the field
  // at the beginning of the simulation. There can be ONLY one initializer.
  std::weak_ptr<FieldInitializer>     m_initializer;

  // If this field is a sub-view of another field, we keep a pointer to the parent
  // On the other hand, if there are sub-views of this field, we keep a list of them
  std::weak_ptr<FieldTracking>              m_parent;
  std::list<std::weak_ptr<FieldTracking>>   m_children;

  // Groups are used to bundle together fields, so that a process can request all of them
  // without knowing/listing all their names. For instance, the dynamics process needs to
  // get all tracers, which need to be advected. However, dyamics has no idea of what are
  // the tracers names, and neither should it care. Groups can come to rescue here, allowing
  // dynamics to request all fields that have been marked as 'tracers'.
  std::set<ci_string>    m_groups;
};

// Use this free function to exploit features of enable_from_this
template<typename... Args>
inline std::shared_ptr<FieldTracking>
create_tracking(const Args&... args) {
  auto ptr = std::make_shared<FieldTracking>(args...);
  ptr->setSelfPointer(ptr);
  return ptr;
}


} // namespace scream

#endif // SCREAM_FIELD_TRACKING_HPP
