#ifndef SCREAM_FIELD_HEADER_HPP
#define SCREAM_FIELD_HEADER_HPP

#include "field_tag.hpp"
#include "util/time_stamp.hpp"

#include <vector>   // For std::vector
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

// A small structure to hold tracking information about a field
// This structure is used to track whether a field has been updated
// (by comparing time stamps), and to track what atm processes
// require/compute a field, so that one can easily check what atm
// process computes a field.
struct FieldTracking {

  // Time stamp of the last update
  util::TimeStamp   m_last_update;

  // List of provider/customer processes. A provider is an atm process
  // that computes/updates the field. A customer is an atm process that
  // uses the field just as an input
  // NOTE: do NOT use shared_ptr, since you will likely create circular references.
  std::vector<std::weak_ptr<AtmosphereProcess>> m_providers;
  std::vector<std::weak_ptr<AtmosphereProcess>> m_customers;
};

// A small class to hold info about a field (everything except the actual field values)
class FieldHeader {
public:

  // Constructor(s)
  FieldHeader () = default;
  FieldHeader (const FieldHeader&) = default;
  FieldHeader (const std::string& name,
               const std::vector<FieldTag>& tags);
  FieldHeader (const std::string& name,
               const std::vector<int>& dims,
               const std::vector<FieldTag>& tags);

  // Assignment (defaulted)
  FieldHeader& operator= (const FieldHeader&) = default;

  // Getters
  const std::string& name () const { return m_name; }
  int rank () const { return m_rank; }
  int dim (const int idim) const {
    error::runtime_check(m_dims.size()==static_cast<size_t>(m_rank), "Error! Field dimensions not yet set.\n",-1);
    error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
    return m_dims[idim];
  } 
  FieldTag tag (const int idim) const { 
    error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
    return m_tags[idim];
  } 
  const FieldTracking& tracking () const { return m_tracking; }

  // If the field was created without specifying the dimensions (e.g., they
  // were not known at construction time), we can reshape the field now.
  // Note: this method cannot be called twice.
  void set_dimensions (const std::vector<int>& dims);

  // The FieldRepository is the only class that can actively modify the tracking
  template<typename MemSpace>
  friend class FieldRepository;

  friend bool operator== (const FieldHeader&, const FieldHeader&);

protected:
  // Some of these could actually be retrieved from the Kokkos View, but it probably makes sens
  // to be able to retrieve info from a header, without having to query the view
  // Besides, we may want to be able to query name/rank/dims BEFORE the view is actually instantiated.
  // Finally, if you look at Field, you may notice that the way we store the view
  // in the manager may lose info about the rank/dims, so we store it here
  std::string           m_name;
  int                   m_rank;
  std::vector<int>      m_dims;
  std::vector<FieldTag> m_tags;

  // Tracking the updates of the field as well as its providers/customers
  FieldTracking         m_tracking;

  // Something about output/restart?
};

bool operator== (const FieldHeader& fh1, const FieldHeader& fh2);
inline bool operator!= (const FieldHeader& fh1, const FieldHeader& fh2) { return !(fh1==fh2); }

} // namespace scream

#endif // SCREAM_FIELD_HEADER_HPP
