#ifndef SCREAM_FIELD_IDENTIFIER_HPP
#define SCREAM_FIELD_IDENTIFIER_HPP

#include "field_tag.hpp"
#include <share/error_defs.hpp>

#include <string>
#include <vector>

namespace scream
{

/*
 *  A small class to hold basic info about a field
 *  
 *  This tiny class is only used to uniquely identify a field, using
 *  the name, the rank, the tags, and possibly the dimensions.
 *  There is no additional meta data about this field.
 */

class FieldIdentifier {
public:

  // Constructor(s)
  FieldIdentifier () = delete;
  FieldIdentifier (const FieldIdentifier&) = default;
  FieldIdentifier (const std::string& name,
                   const std::vector<FieldTag>& tags);

  // Assignment (defaulted)
  FieldIdentifier& operator= (const FieldIdentifier&) = default;

  // ----- Getters ----- //

  // Name and layout informations
  const std::string& name () const { return m_name; }
  const std::vector<FieldTag>& tags () const { return m_tags; }
  FieldTag tag  (const int idim) const;
  int      rank ()               const  { return m_rank; }
  const std::vector<int>& dims () const { return m_dims; }
  int      dim  (const int idim) const;
  int      size ()               const;

  // The identifier string
  const std::string& get_identifier () const { return m_identifier; }

  bool is_dimension_set  (const int idim) const;
  bool are_dimensions_set () const;
  
  // ----- Setters ----- //

  // Note: as soon as a dimension is set, it cannot be changed.
  void set_dimension  (const int idim, const int dimension);
  void set_dimensions (const std::vector<int>& dims);

  // We reimplement the equality operator for identifiers comparison (needed for some std container)
  friend bool operator== (const FieldIdentifier&, const FieldIdentifier&);
  friend bool operator<  (const FieldIdentifier&, const FieldIdentifier&);

protected:

  void update_identifier ();

  std::string           m_name;
  std::vector<FieldTag> m_tags;
  int                   m_rank;
  std::vector<int>      m_dims;

  // The identifier string is a conveniet way to display the information of
  // the identifier, so that it can be easily read.
  std::string           m_identifier;
};

bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2);
inline bool operator!= (const FieldIdentifier& fid1, const FieldIdentifier& fid2) { return !(fid1==fid2); }

// ========================== IMPLEMENTATION ======================= //

inline int FieldIdentifier::dim (const int idim) const {
  error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_dims[idim];
}

inline int FieldIdentifier::size () const {
  error::runtime_check(are_dimensions_set(), "Error! Field dimensions not yet set.\n",-1);
  int prod = m_rank>0 ? 1 : 0;
  for (int idim=0; idim<m_rank; ++idim) {
    prod *= m_dims[idim];
  }
  return prod;
}

inline FieldTag FieldIdentifier::tag (const int idim) const { 
  error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_tags[idim];
} 

inline bool FieldIdentifier::is_dimension_set (const int idim) const {
  error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_dims[idim]>=0;
}

inline bool FieldIdentifier::are_dimensions_set () const {
  for (int idim=0; idim<m_rank; ++idim) {
    if (m_dims[idim]<0) {
      return false;
    }
  }
  return true;
}

} // namespace scream

#endif // SCREAM_FIELD_IDENTIFIER_HPP
