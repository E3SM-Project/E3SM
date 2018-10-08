#ifndef SCREAM_FIELD_HEADER_HPP
#define SCREAM_FIELD_HEADER_HPP

#include "field_tag.hpp"
#include "util/time_stamp.hpp"

#include <set>      // For std::set
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

// A small class to hold info about a field (everything except the actual field values)
class FieldHeader {
public:

  // The type of the customer/provider sets.
  // NOTE: weak_ptr does not implement an operator<, so we need this long name
  using atm_procs_set = std::set<std::weak_ptr<AtmosphereProcess>,std::owner_less<std::weak_ptr<AtmosphereProcess>>>;

  // Constructor(s)
  FieldHeader () = default;
  FieldHeader (const FieldHeader&) = default;
  FieldHeader (const std::string& name,
               const std::vector<FieldTag>& tags);
  FieldHeader (const std::string& name,
               const std::vector<FieldTag>& tags,
               const std::vector<int>& dims);

  // Assignment (defaulted)
  FieldHeader& operator= (const FieldHeader&) = default;

  // Getters
  const std::string& name () const { return m_name; }
  FieldTag tag  (const int idim) const;
  int      dim  (const int idim) const;
  int      size ()               const;
  int      rank ()               const  { return m_rank; }
  // The time stamp of the field. This can be used to check when it was last update.
  // Please, notice this is not the OS time stamp (see TimeStamp.hpp for details).
  const util::TimeStamp& time_stamp    () const { return m_time_stamp; }
  const std::string&     identifier    () const { return m_identifier; }
  const atm_procs_set&   get_providers () const { return m_providers; }
  const atm_procs_set&   get_customers () const { return m_customers; }

  // Setters
  // If the field was created without specifying the dimensions (e.g., they
  // were not known at construction time), we can reshape the field now.
  // Note: this method cannot be called twice.
  void set_dimensions (const std::vector<int>& dims);
  bool dimensions_set () const { return m_dims.size()==static_cast<size_t>(m_rank); }

  // Add to the list of providers/customers
  void add_provider (std::shared_ptr<AtmosphereProcess> provider);
  void add_customer (std::shared_ptr<AtmosphereProcess> customer);

  friend bool operator== (const FieldHeader&, const FieldHeader&);

protected:

  void update_identifier ();

  // Some of these could actually be retrieved from the Kokkos View, but it probably makes sens
  // to be able to retrieve info from a header, without having to query the view
  // Besides, we may want to be able to query name/rank/dims BEFORE the view is actually instantiated.
  // Finally, if you look at Field, you may notice that the way we store the view
  // in the manager may lose info about the rank/dims, so we store it here
  std::string           m_name;
  std::vector<FieldTag> m_tags;
  std::vector<int>      m_dims;
  int                   m_rank;

  // A string used to identify this field. This is more than just the field name,
  // since we can have different fields with the same name (e.g., with different layouts).
  // The identifier puts the name, tags and dims info in one string.
  // The only purpose of this string is to put all fields into one std::map.
  std::string           m_identifier;

  // Tracking the updates of the field
  util::TimeStamp   m_time_stamp;

  // List of provider/customer processes. A provider is an atm process
  // that computes/updates the field. A customer is an atm process that
  // uses the field just as an input
  // NOTE: do NOT use shared_ptr, since you will likely create circular references.
  atm_procs_set m_providers;
  atm_procs_set m_customers;
};

bool operator== (const FieldHeader& fh1, const FieldHeader& fh2);
inline bool operator!= (const FieldHeader& fh1, const FieldHeader& fh2) { return !(fh1==fh2); }

// ========================== IMPLEMENTATION ======================= //

inline FieldTag FieldHeader::tag (const int idim) const { 
  error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_tags[idim];
} 

inline int FieldHeader::dim (const int idim) const {
  error::runtime_check(m_dims.size()==static_cast<size_t>(m_rank), "Error! Field dimensions not yet set.\n",-1);
  error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_dims[idim];
}

inline int FieldHeader::size () const {
  error::runtime_check(m_dims.size()==static_cast<size_t>(m_rank), "Error! Field dimensions not yet set.\n",-1);
  int prod = m_rank>0 ? 1 : 0;
  for (int idim=0; idim<m_rank; ++idim) {
    prod *= m_dims[idim];
  }
  return prod;
}

} // namespace scream

#endif // SCREAM_FIELD_HEADER_HPP
