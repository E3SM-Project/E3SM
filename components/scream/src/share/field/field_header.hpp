#ifndef SCREAM_FIELD_HEADER_HPP
#define SCREAM_FIELD_HEADER_HPP

#include "field_identifier.hpp"
#include "field_tracking.hpp"
#include "field_tag.hpp"
#include <share/scream_types.hpp>
#include <share/util/time_stamp.hpp>

#include <vector>
#include <memory>   // For std::shared_ptr and std::weak_ptr

namespace scream
{

class AtmosphereProcess;

// A small class to hold info about a field (everything except the actual field values)
class FieldHeader {
public:

  using identifier_type = FieldIdentifier;
  using tracking_type   = FieldTracking;

  // Constructor(s)
  FieldHeader () = default;
  FieldHeader (const FieldHeader&) = default;
  FieldHeader (const identifier_type& id);

  // Assignment (defaulted)
  FieldHeader& operator= (const FieldHeader&) = default;

  // ----- Getters ----- //

  // Get the basic information from the identifier
  const FieldIdentifier& get_identifier () const { return m_identifier; }

  // Get the tracking
  const FieldTracking& get_tracking () const { return m_tracking; }
        FieldTracking& get_tracking ()       { return m_tracking; }

protected:

  // Static information about the field: name, rank, tags
  FieldIdentifier       m_identifier;

  // Tracking of the field
  FieldTracking         m_tracking;
};

} // namespace scream

#endif // SCREAM_FIELD_HEADER_HPP
