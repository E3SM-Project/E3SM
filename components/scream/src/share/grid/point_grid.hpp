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
 * This grid is typical of Physics parametrizations, and is also used
 * to interface to the component coupler.
 */

class PointGrid : public AbstractGrid
{
public:

  PointGrid (const std::string& grid_name,
             const int num_global_cols,
             const int num_my_cols,
             const int num_vertical_levels);

  virtual ~PointGrid () = default;

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a 2d structured grid, this could be a set of 2 indices.
  FieldLayout get_native_dof_layout () const override;

  void set_dofs (const dofs_list_type& dofs);
  void set_geometry_data (const std::string& name, const geo_view_type& data) override;
};

// Create a point grid, with linear range of gids, evenly partitioned
// among the ranks in the given communicator.
std::shared_ptr<const PointGrid>
create_point_grid (const std::string& name,
                   const int num_global_cols,
                   const int num_vertical_lev,
                   const ekat::Comm& comm);

} // namespace scream

#endif // SCREAM_POINT_GRID_HPP
