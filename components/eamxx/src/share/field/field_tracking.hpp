#ifndef SCREAM_FIELD_TRACKING_HPP
#define SCREAM_FIELD_TRACKING_HPP

#include "share/field/field_group_info.hpp"
#include "share/eamxx_types.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/util/eamxx_family_tracking.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <memory>   // For std::weak_ptr
#include <string>

namespace scream {

// Forward declarations
class AtmosphereProcess;
class FieldHeader;

class FieldTracking : public FamilyTracking<FieldTracking> {
public:

  using TimeStamp         = util::TimeStamp;
  using ci_string         = ekat::CaseInsensitiveString;
  using atm_proc_ptr_type = std::weak_ptr<AtmosphereProcess>;
  using atm_proc_set_type = ekat::WeakPtrSet<AtmosphereProcess>;

  FieldTracking() = default;
  FieldTracking(const FieldTracking&) = default;

  // No assignment, to prevent tampering with tracking (e.g., rewinding time stamps)
  FieldTracking& operator=(const FieldTracking&) = delete;

  // ----- Getters ----- //

  // The time stamp of the field. This can be used to check when it was last updated.
  // Please, notice this is not the OS time stamp (see time_stamp.hpp for details).
  const TimeStamp& get_time_stamp () const { return m_time_stamp; }

  //  - provider: can compute the field as an output
  //  - customer: requires the field as an input
  const atm_proc_set_type& get_providers () const { return m_providers; }
  const atm_proc_set_type& get_customers () const { return m_customers; }

  // List of field groups that this field belongs to
  const ekat::WeakPtrSet<const FieldGroupInfo>& get_groups_info () const { return m_groups; }

  // ----- Setters ----- //

  // Add to the list of providers/customers
  void add_provider (const std::weak_ptr<AtmosphereProcess>& provider);
  void add_customer (const std::weak_ptr<AtmosphereProcess>& customer);

  // Add the field to a given group
  void add_to_group (const std::shared_ptr<const FieldGroupInfo>& group);

  // Set the time stamp for this field. This can only be called once, due to TimeStamp implementation.
  // NOTE: if the field has 'children' (see FamilyTracking), their ts will be updated too.
  //       However, if the field has a 'parent' (see FamilyTracking), the parent's ts will not be updated.
  void update_time_stamp (const TimeStamp& ts);
  void invalidate_time_stamp ();

  // Set/get accumulation interval start
  void set_accum_start_time (const TimeStamp& ts);
  const TimeStamp& get_accum_start_time () const { return m_accum_start; }

protected:

  // We keep the field name just to make debugging messages more helpful
  std::string m_name;

  // Tracking the updates of the field
  TimeStamp         m_time_stamp;

  // For accumulated vars, the time where the accumulation started
  TimeStamp         m_accum_start;
  ci_string         m_accum_type;

  // List of provider/customer processes. A provider is an atm process that computes/updates the field.
  // A customer is an atm process that uses the field just as an input.
  // NOTE: do NOT use shared_ptr, since you would create circular references.
  atm_proc_set_type       m_providers;
  atm_proc_set_type       m_customers;

  // Groups are used to bundle together fields, so that a process can request all of them
  // without knowing/listing all their names. For instance, the dynamics process needs to
  // get all tracers, which need to be advected. However, dynamics has no idea (a priori)
  // of what are the tracers names or how many there are, and neither should it care.
  // FieldGroup's allow atm procs to request all fields that have been marked as 'tracers'.
  // Here, we keep track of all the groups that this field belongs to.
  ekat::WeakPtrSet<const FieldGroupInfo>    m_groups;
};

// Use this free function to exploit features of enable_shared_from_this,
// as well as features from FamilyTracking.
template<typename... Args>
inline std::shared_ptr<FieldTracking>
create_tracking(const Args&... args) {
  auto ptr = std::make_shared<FieldTracking>(args...);
  ptr->setSelfPointer(ptr);
  return ptr;
}


} // namespace scream

#endif // SCREAM_FIELD_TRACKING_HPP
