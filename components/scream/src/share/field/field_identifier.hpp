#ifndef SCREAM_FIELD_IDENTIFIER_HPP
#define SCREAM_FIELD_IDENTIFIER_HPP

#include "share/field/field_layout.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/scream_assert.hpp"

#include <string>
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
  using layout_type   = FieldLayout;
  using grid_type     = AbstractGrid;
  using grid_ptr_type = std::shared_ptr<grid_type>;

  // Constructor(s)
  FieldIdentifier () = delete;
  FieldIdentifier (const FieldIdentifier&) = default;
  FieldIdentifier (const std::string& name,
                   const layout_type& layout,
                   const grid_ptr_type grid = nullptr);
  FieldIdentifier (const std::string& name,
                   const std::vector<FieldTag>& tags,
                   const grid_ptr_type grid = nullptr);
  FieldIdentifier (const std::string& name,
                   const std::initializer_list<FieldTag>& tags,
                   const grid_ptr_type grid = nullptr);

  // Assignment (defaulted)
  FieldIdentifier& operator= (const FieldIdentifier&) = default;

  // ----- Getters ----- //

  // Name and layout informations
  const std::string&  name       () const { return m_name;   }
  const layout_type&  get_layout () const { return m_layout; }
        grid_ptr_type get_grid   () const { return m_grid;   }

  // The identifier string
  const std::string& get_identifier () const { return m_identifier; }

  // ----- Setters ----- //

  // Note: as soon as a dimension is set, it cannot be changed.
  void set_dimension  (const int idim, const int dimension);
  void set_dimensions (const std::vector<int>& dims);
  void set_grid       (const grid_ptr_type grid);

  // We reimplement the equality operator for identifiers comparison (needed for some std container)
  friend bool operator== (const FieldIdentifier&, const FieldIdentifier&);
  friend bool operator<  (const FieldIdentifier&, const FieldIdentifier&);

protected:

  void update_identifier ();

  std::string     m_name;

  layout_type     m_layout;

  grid_ptr_type   m_grid;

  // The identifier string is a conveniet way to display the information of
  // the identifier, so that it can be easily read.
  std::string     m_identifier;
};

bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2);
inline bool operator!= (const FieldIdentifier& fid1, const FieldIdentifier& fid2) { return !(fid1==fid2); }

} // namespace scream

#endif // SCREAM_FIELD_IDENTIFIER_HPP
