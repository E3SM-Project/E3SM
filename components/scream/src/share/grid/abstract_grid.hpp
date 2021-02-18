#ifndef SCREAM_ABSTRACT_GRID_HPP
#define SCREAM_ABSTRACT_GRID_HPP

#include "share/grid/grid_utils.hpp"
#include "share/field/field_layout.hpp"
#include "share/scream_types.hpp"

#include <map>

namespace scream
{

/*
 * An interface base class for Grid objects
 *
 * An abstract grid represents a 2d grid, vertically extruded. It needs to be able
 * to provide the following information:
 *   - number of local degrees of freedom (dofs) along the horizontal direction
 *   - number of vertical levels
 *   - gids for all local (2d) dofs
 *   - type and name of the grid (see grid_utils.hpp for types)
 *   - the layout of a (2d) dof on this grid (see field_layout.hpp)
 *   - the gid of a (2d) dof given indices (i_1,...,i_N), with N being the
 *     rank of the scalar field layout, and i_k less than the i-th dimension
 *     in the scalar field layout
 *
 * The native dof layout is the way one would naturally index dofs in the grid.
 * For instance, on a cartesian 2d grid, this would be done with 2 indices.
 * On a Spectral Element grid, it would be 3 indices (element, gauss pt, gauss pt).
 * This is the layout that a 2d scalar field would have on this grid.
 * Having the grid exposing this layout, allows downstream classes to not avoid
 * assumptions or switches on the particular grid, and simply rely on common
 * grid interfaces.
 */

class AbstractGrid
{
public:
  using gid_type          = int;           // TODO: use int64_t? int? template class on gid_type?
  using device_type       = DefaultDevice; // TODO: template class on device type
  using kokkos_types      = KokkosTypes<device_type>;
  using geo_view_type     = kokkos_types::view_1d<double>;
  using geo_view_map_type = std::map<std::string,geo_view_type>;

  // The list of all dofs' gids
  using dofs_list_type = kokkos_types::view_1d<gid_type>;

  // Row i of this 2d view gives the indices of the ith local dof
  // in the native 2d layout
  using lid_to_idx_map_type = kokkos_types::view<int**>;

  // Constructor(s) & Destructor
  AbstractGrid (const GridType type,
                const std::string& name)
   : AbstractGrid(0,0,type,name)
  {
    // Nothing to do here
  }

  AbstractGrid (int num_local_dofs,
                int num_global_dofs,
                const GridType type,
                const std::string& name)
   : m_type (type)
   , m_name (name)
   , m_num_local_dofs  (num_local_dofs)
   , m_num_global_dofs (num_global_dofs)
  {
    // Nothing to do here
  }

  virtual ~AbstractGrid () = default;

  // Grid description utilities
  GridType type () const { return m_type; }
  const std::string& name () const { return m_name; }

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a scalar 2d field on a SE grid, this will be (nelem,np,np),
  //       for a vector 3d field on a Point grid it will be (ncols,vector_dim,nlevs)
  virtual FieldLayout get_2d_scalar_layout () const = 0;
  virtual FieldLayout get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const = 0;
  virtual FieldLayout get_3d_scalar_layout (const bool midpoints) const = 0;
  virtual FieldLayout get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const = 0;

  int get_num_vertical_levels () const { return m_num_vert_levs; }

  // The number of dofs on this MPI rank
  int get_num_local_dofs  () const { return m_num_local_dofs;  }
  gid_type get_num_global_dofs () const { return m_num_global_dofs; }

  // Get a 1d view containing the dof gids
  const dofs_list_type& get_dofs_gids () const { return m_dofs_gids; }

  // Get a 2d view, where (i,j) entry contains the j-th coordinate of
  // the i-th dof in the native dof layout.
  const lid_to_idx_map_type& get_lid_to_idx_map () const { return m_lid_to_idx; }

  // Set/get geometric views. The setter is virtual, so each grid can check if "name" is supported.
  virtual void set_geometry_data (const std::string& name, const geo_view_type& data) = 0;
  const geo_view_type& get_geometry_data (const std::string& name) const;

protected:

  // The grid name and type
  GridType     m_type;
  std::string  m_name;

  // Counters
  int m_num_local_dofs;
  int m_num_global_dofs;
  int m_num_vert_levs;

  // The global ID of each dof
  dofs_list_type        m_dofs_gids;

  // The map lid->idx
  lid_to_idx_map_type   m_lid_to_idx;

  geo_view_map_type     m_geo_views;
};

inline const AbstractGrid::geo_view_type&
AbstractGrid::get_geometry_data (const std::string& name) const {
  EKAT_REQUIRE_MSG (m_geo_views.find(name)!=m_geo_views.end(),
                    "Error! Grid '" + m_name + "' does not store geometric data '" + name + "'.\n");
  return m_geo_views.at(name);
}

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
