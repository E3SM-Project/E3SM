#ifndef EAMXX_LATLON_GRID_HPP
#define EAMXX_LATLON_GRID_HPP

#include "ekat/mpi/ekat_comm.hpp"
#include "share/eamxx_types.hpp"
#include "share/grid/abstract_grid.hpp"

namespace scream {

/*
 * A grid consisting of a bunch of points structured on a lat-lon grid.
 *
 * A generalization of the point grid impl.
 *
 */

class LatLonGrid : public AbstractGrid {
 public:
  LatLonGrid(const std::string &grid_name, const int nglat, const int nglon,
             const int my_nlon, const int num_vertical_levels,
             const ekat::Comm &comm);

  virtual ~LatLonGrid() = default;

  // Native layout of a dof. This is the natural way to index a dof in the grid,
  // e.g., for a 2d structured grid, this could be a set of 2 indices.
  FieldLayout get_2d_scalar_layout() const override;
  FieldLayout get_2d_vector_layout(
      const int vector_dim, const std::string &vec_dim_name) const override;
  FieldLayout get_2d_tensor_layout(
      const std::vector<int> &cmp_dims,
      const std::vector<std::string> &cmp_names) const override;
  FieldLayout get_3d_scalar_layout(const bool midpoints) const override;
  FieldLayout get_3d_vector_layout(
      const bool midpoints, const int vector_dim,
      const std::string &vec_dim_name) const override;
  FieldLayout get_3d_tensor_layout(
      const bool midpoints, const std::vector<int> &cmp_dims,
      const std::vector<std::string> &cmp_names) const override;

  FieldTag get_partitioned_dim_tag() const override {
    return FieldTag::Longitude;
  }
  int get_partitioned_dim_local_size() const override { return m_nlon_local; }
  int get_partitioned_dim_global_size() const override { return m_nlon_global; }

  std::shared_ptr<AbstractGrid> clone(const std::string &clone_name,
                                      const bool shallow) const override;

 protected:
  int m_nlat_global;
  int m_nlon_global;
  int m_nlon_local;
};

// Create a lat-lon grid
// with gids evenly partitioned among the ranks in the given communicator.
std::shared_ptr<LatLonGrid> create_latlon_grid(const std::string &name,
                                               const int nglat, const int nglon,
                                               const int num_vertical_lev,
                                               const ekat::Comm &comm);

}  // namespace scream

#endif  // EAMXX_LATLON_GRID_HPP
