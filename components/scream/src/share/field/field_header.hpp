#ifndef SCREAM_FIELD_HEADER_HPP
#define SCREAM_FIELD_HEADER_HPP

#include "field_identifier.hpp"
#include "field_tag.hpp"
#include <share/util/time_stamp.hpp>

#include <vector>
#include <memory>   // For std::shared_ptr and std::weak_ptr

// I have to decide where to store a field's providers and customers, that is,
// a list of parametrizations that compute or need a given field. Should it
// be done by the FieldRepository? Or should a Field store its own? I lean toward
// the latter.
// NOTE: if two parametrizations both compute/update a field, we need to ensure 
//       that the two parametrizations are run sequentially. The FieldTracking
//       structure can be used (among other things) to enforce this requirement,
//       by checking the time stamps of the fields.

namespace scream
{

class AtmosphereProcess;

// A small class to hold info about a field (everything except the actual field values)
class FieldHeader {
public:

  using identifier_type = FieldIdentifier;

  // Constructor(s)
  FieldHeader () = default;
  FieldHeader (const FieldHeader&) = default;
  FieldHeader (const identifier_type& id);

  // Assignment (defaulted)
  FieldHeader& operator= (const FieldHeader&) = default;

  // ----- Getters ----- //

  // Get the basic information from the identifier
  const FieldIdentifier& get_identifier () const { return m_identifier; }

  // The time stamp of the field. This can be used to check when it was last update.
  // Please, notice this is not the OS time stamp (see TimeStamp.hpp for details).
  const util::TimeStamp& get_time_stamp    () const { return m_time_stamp; }

  // List of providers/customers for this field
  const std::vector<std::shared_ptr<AtmosphereProcess>>& get_providers () const { return m_providers; }
  const std::vector<std::shared_ptr<AtmosphereProcess>>& get_customers () const { return m_customers; }

  // List of field groups that this field belongs to
  const std::vector<std::string>& get_groups_list () const { return m_groups; }

  // ----- Setters ----- //

  // If the header was created without specifying dimensions
  // (e.g., they were not known at construction time),
  // call this method to set them.
  void set_dimensions (const std::vector<int>& dims) { m_identifier.set_dimensions(dims); }

  // Add to the list of providers/customers
  void add_provider (std::shared_ptr<AtmosphereProcess> provider);
  void add_customer (std::shared_ptr<AtmosphereProcess> customer);

  // Add the field to a given group
  void add_to_group (const std::string& group_name);

  // We reimplement the equality operator for headers comparison (needed for std::set)
  // friend bool operator== (const FieldHeader&, const FieldHeader&);

protected:

  // Static information about the field: name, rank, and tags
  FieldIdentifier           m_identifier;

  // Tracking the updates of the field
  util::TimeStamp   m_time_stamp;

  // These are to be used to track the order in which providers update the field at each time step.
  // One can use this information to track when a field gets updated during a timestep. It can be
  // particularly useful in the case of parallel schedules.
  std::vector<std::string> m_last_ts_providers;
  std::vector<std::string> m_curr_ts_providers;

  // List of provider/customer processes. A provider is an atm process
  // that computes/updates the field. A customer is an atm process that
  // uses the field just as an input
  // NOTE: do NOT use shared_ptr, since you will likely create circular references.
  std::vector<std::shared_ptr<AtmosphereProcess>> m_providers;
  std::vector<std::shared_ptr<AtmosphereProcess>> m_customers;

  // Groups are used to bundle together fields, so that a process can request all of them
  // without knowing/listing all their names. For instance, the dynamics process needs to
  // get all tracers, which need to be advected. However, dyamics has no idea of what are
  // the tracers names, and neither should it care. Groups can come to rescue here, allowing
  // dynamics to request all fields that have been marked as 'tracers'.
  std::vector<std::string>    m_groups;
};

} // namespace scream

#endif // SCREAM_FIELD_HEADER_HPP
