#include "share/field/field_identifier.hpp"
#include "share/grid/default_grid.hpp"

namespace scream
{

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const grid_ptr_type grid)
 : m_name   (name)
 , m_layout (layout)
{
  // This also calls 'update_identifier'
  set_grid(grid);
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const std::vector<FieldTag>& tags,
                 const grid_ptr_type grid)
 : m_name   (name)
 , m_layout (tags)
{
  // This also calls 'update_identifier'
  set_grid(grid);
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const std::initializer_list<FieldTag>& tags,
                 const grid_ptr_type grid)
 : m_name   (name)
 , m_layout (tags)
{
  // This also calls 'update_identifier'
  set_grid(grid);
}

void FieldIdentifier::set_dimension (const int idim, const int dimension) {
  m_layout.set_dimension(idim,dimension);
  update_identifier ();
}

void FieldIdentifier::set_dimensions (const std::vector<int>& dims) {
  m_layout.set_dimensions(dims);
  update_identifier ();
}

void FieldIdentifier::set_grid (const grid_ptr_type grid) {
  // Only allow overwriting if the stored grid is empty
  scream_require_msg (!static_cast<bool>(m_grid) || m_grid->type()==GridType::Undefined, "Error! Cannot overwrite a non-empty grid.\n");
  if (grid==nullptr) {
    // create empty grid and set that
    m_grid = std::make_shared<DefaultGrid<GridType::Undefined>>();
  } else {
    m_grid = grid;
  }

  // Update the identifier string
  update_identifier ();
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier = m_name + "[" + e2str(m_grid->type()) + "]<" + tag2string(m_layout.tags()[0]);
  for (int dim=1; dim<m_layout.rank(); ++dim) {
    m_identifier += "," + tag2string(m_layout.tags()[dim]);
  }
  m_identifier += ">(" + std::to_string(m_layout.dims()[0]);
  for (int dim=1; dim<m_layout.rank(); ++dim) {
    m_identifier += "," + std::to_string(m_layout.dims()[dim]);
  }
  m_identifier += ")";
}

// Free functions for identifiers comparison
bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_identifier()==fid2.get_identifier());
}

bool operator< (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_identifier()<fid2.get_identifier());
}

} // namespace scream
