#ifndef SCREAM_POINT_GRID_HPP
#define SCREAM_POINT_GRID_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"

#include "ekat/mpi/ekat_comm.hpp"

namespace scream
{

/*
 * A grid consisting of a bunch of points
 *
 * A point grid simply stores a bunch of dofs gids. Unlike SEGrid,
 * which also stores info on the SE element each dof belongs to,
 * this class does not store any additional info regarding its dofs.
 * In particular, the map lid->idx is an identity, and the native
 * layout has only one tag: Column.
 *
 * This grid is typical of Physics parametrizations
 */

class PointGrid : public AbstractGrid
{
public:

  PointGrid (const std::string& grid_name,
             const dofs_list_type& dofs_gids,
             const int num_vertical_lev);
  virtual ~PointGrid () = default;

  // Grid description utilities
  GridType type () const override { return GridType::Point; }
  const std::string& name () const override { return m_grid_name; }

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a 2d structured grid, this could be a set of 2 indices.
  FieldLayout get_native_dof_layout () const override;

  int get_num_vertical_levels () const override { return m_num_vl; }

  // Dofs gids utilities
  int get_num_local_dofs () const override { return m_num_my_cols; }
  const dofs_list_type& get_dofs_gids () const override { return m_dofs_gids; }
  lid_to_idx_map_type get_lid_to_idx_map () const override { return m_lid_to_idx; }

protected:

  const std::string     m_grid_name;

  int                   m_num_my_cols;
  int                   m_num_vl;
  dofs_list_type        m_dofs_gids;
  lid_to_idx_map_type   m_lid_to_idx;
};

// Create a point grid, with linear range of gids, evenly partitioned
// among the ranks in the given communicator.
PointGrid
create_point_grid (const std::string& name,
                   const int num_global_cols,
                   const int num_vertical_lev,
                   const ekat::Comm& comm);

} // namespace scream

#endif // SCREAM_POINT_GRID_HPP
