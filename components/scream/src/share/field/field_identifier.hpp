#ifndef SCREAM_FIELD_IDENTIFIER_HPP
#define SCREAM_FIELD_IDENTIFIER_HPP

#include "share/field/field_layout.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/util/ekat_units.hpp"

#include <vector>

namespace scream
{

/*
 *  A small class to hold basic info about a field
 *
 *  This tiny class is only used to uniquely identify a field, using
 *  the name, and its layout (which containss the rank, the tags,
 *  and possibly the dimensions).
 *  There is no additional meta data about this field.
 */

class FieldIdentifier {
public:
  using layout_type = FieldLayout;
  using ci_string   = ekat::CaseInsensitiveString;
  using Units       = ekat::units::Units;

  // Constructor(s)
  FieldIdentifier () = delete;
  FieldIdentifier (const FieldIdentifier&) = default;
  FieldIdentifier (const std::string& name,
                   const layout_type& layout,
                   const Units& units,
                   const std::string& grid_name = "");
  FieldIdentifier (const std::string& name,
                   const std::vector<FieldTag>& tags,
                   const Units& units,
                   const std::string& grid_name = "");
  FieldIdentifier (const std::string& name,
                   const std::initializer_list<FieldTag>& tags,
                   const Units& units,
                   const std::string& grid_name = "");

  // Delete assignment, to prevent overwriting identifiers sneakyly
  FieldIdentifier& operator= (const FieldIdentifier&) = delete;

  // ----- Getters ----- //

  // Name and layout informations
  const std::string&  name          () const { return m_name;      }
  const layout_type&  get_layout    () const { return m_layout;    }
  const Units&        get_units     () const { return m_units;     }
  const std::string&  get_grid_name () const { return m_grid_name; }

  // The identifier string
  const std::string& get_id_string () const { return m_identifier; }

  // ----- Setters ----- //

  // Note: as soon as a dimension is set, it cannot be changed.
  void set_dimension  (const int idim, const int dimension);
  void set_dimensions (const std::vector<int>& dims);
  void set_grid_name  (const std::string& grid_name);

  // We reimplement the equality operator for identifiers comparison (needed for some std container)
  friend bool operator== (const FieldIdentifier&, const FieldIdentifier&);
  friend bool operator<  (const FieldIdentifier&, const FieldIdentifier&);

protected:

  void update_identifier ();

  ci_string       m_name;

  layout_type     m_layout;

  Units           m_units;

  ci_string       m_grid_name;

  // The identifier string is a conveniet way to display the information of
  // the identifier, so that it can be easily read.
  ci_string       m_identifier;
};

bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2);
inline bool operator!= (const FieldIdentifier& fid1, const FieldIdentifier& fid2) { return !(fid1==fid2); }

} // namespace scream

#endif // SCREAM_FIELD_IDENTIFIER_HPP
