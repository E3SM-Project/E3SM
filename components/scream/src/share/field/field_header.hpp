#ifndef SCREAM_FIELD_HEADER_HPP
#define SCREAM_FIELD_HEADER_HPP

#include "share/field/field_identifier.hpp"
#include "share/field/field_tracking.hpp"
#include "share/field/field_alloc_prop.hpp"
#include "share/scream_types.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "ekat/std_meta/ekat_std_any.hpp"

#include <vector>
#include <map>
#include <memory>   // For std::shared_ptr and std::weak_ptr

namespace scream
{

class AtmosphereProcess;

/*
 * A small class to contain meta-data about a field
 *
 * The FieldHeader class is itself a container of other
 * more speicific classes, such as FieldIdentifier
 * (which contains information used to uniquely identify
 * the field) or FieldTracking (which contains info used
 * to track access to the field).
 * There is also 'extra_data', which is a sort of fall-back
 * option, for the meta-data that does not follow under
 * any pre-defined category, and that is not general enough
 * to warrant a new sub-object or a specific named member/method.
 */

class FieldHeader {
public:

  using identifier_type = FieldIdentifier;
  using tracking_type   = FieldTracking;
  using extra_data_type = std::map<std::string,ekat::util::any>;

  // Constructor(s)
  FieldHeader (const FieldHeader&) = default;
  explicit FieldHeader (const identifier_type& id);

  // Assignment deleted, to prevent sneaky overwrites.
  FieldHeader& operator= (const FieldHeader&) = delete;

  // Set extra data
  void set_extra_data (const std::string& key,
                       const ekat::util::any& data,
                       const bool throw_if_existing = false);

  template<typename T>
  void set_extra_data (const std::string& key,
                       const T& data,
                       const bool throw_if_existing = false) {
    ekat::util::any data_any;
    data_any.reset<T>(data);
    set_extra_data(key,data_any,throw_if_existing);
  }

  // ----- Getters ----- //

  // Get the basic information from the identifier
  const FieldIdentifier& get_identifier () const { return m_identifier; }

  // Get the tracking
  const FieldTracking& get_tracking () const { return m_tracking; }
        FieldTracking& get_tracking ()       { return m_tracking; }

  // Get the allocation properties
  const FieldAllocProp& get_alloc_properties () const { return m_alloc_prop; }
        FieldAllocProp& get_alloc_properties ()       { return m_alloc_prop; }

  // Get the extra data
  const extra_data_type& get_extra_data () const { return m_extra_data; }

protected:

  // Static information about the field: name, rank, tags
  FieldIdentifier       m_identifier;

  // Tracking of the field
  FieldTracking         m_tracking;

  // Allocation properties
  FieldAllocProp        m_alloc_prop;

  // Extra data associated with this field
  extra_data_type       m_extra_data;
};

} // namespace scream

#endif // SCREAM_FIELD_HEADER_HPP
