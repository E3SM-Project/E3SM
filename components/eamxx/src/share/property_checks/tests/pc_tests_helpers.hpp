#include "share/grid/point_grid.hpp"

#include <ekat_units.hpp>

namespace scream {

inline std::shared_ptr<const AbstractGrid>
create_test_grid (const ekat::Comm& comm, int num_lcols, int nlevs)
{
  using namespace ekat::units;
  using gid_type = AbstractGrid::gid_type;

  // Create a point grid
  const auto grid = create_point_grid("some_grid",num_lcols*comm.size(),nlevs,comm);
  const auto layout = grid->get_2d_scalar_layout();
  const auto units = ekat::units::Units::nondimensional();
  const auto& lat = grid->create_geometry_data("lat",layout,units);
  const auto& lon = grid->create_geometry_data("lon",layout,units);
  auto lat_h = lat.get_strided_view<Real*,Host>();
  auto lon_h = lon.get_strided_view<Real*,Host>();
  auto dofs = grid->get_dofs_gids();
  auto dofs_h = dofs.get_strided_view<gid_type*,Host>();
  for (int i=0; i<grid->get_num_local_dofs(); ++i) {
    lat_h(i) = i;
    lon_h(i) = -i;
  }
  lat.sync_to_dev();
  lon.sync_to_dev();

  return grid;
}

inline Field create_test_field (const std::shared_ptr<const AbstractGrid>& grid)
{
  using namespace ekat::units;
  const auto layout = grid->get_3d_vector_layout(true,3);

  FieldIdentifier fid ("field_1", layout, m/s, grid->name());
  Field f(fid);
  f.allocate_view();
  return f;
}

inline Field create_data_field (const std::shared_ptr<const AbstractGrid>& grid)
{
  using namespace ekat::units;
  const auto layout = grid->get_2d_scalar_layout();

  FieldIdentifier fid ("data", layout, m/s, grid->name());
  Field f(fid);
  f.allocate_view();
  f.deep_copy(1);
  return f;
}

} // namespace scream
