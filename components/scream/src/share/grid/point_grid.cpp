#include "share/grid/point_grid.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/share/physics_constants.hpp"

#include <numeric>

namespace scream {

PointGrid::
PointGrid (const std::string& grid_name,
           const int          num_my_cols,
           const int          num_vertical_levels,
           const ekat::Comm&  comm)
 : AbstractGrid(grid_name,GridType::Point,num_my_cols,num_vertical_levels,comm)
{
  // The lid->idx map is the identity map.
  lid_to_idx_map_type lid_to_idx("lid to idx",get_num_local_dofs(),1);
  auto h_lid_to_idx = Kokkos::create_mirror_view(lid_to_idx);
  std::iota(h_lid_to_idx.data(),h_lid_to_idx.data()+get_num_local_dofs(),0);
  Kokkos::deep_copy(lid_to_idx,h_lid_to_idx);

  // Note: we use the base class set method, rather than directly using m_lid_to_idx,
  //       so that we can perform sanity checks on inputs consistency.
  set_lid_to_idx_map(lid_to_idx);
}


FieldLayout
PointGrid::get_2d_scalar_layout () const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({COL},{get_num_local_dofs()});
}

FieldLayout
PointGrid::get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({COL,vector_tag},{get_num_local_dofs(),vector_dim});
}

FieldLayout
PointGrid::get_3d_scalar_layout (const bool midpoints) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({COL,VL},{get_num_local_dofs(),nvl});
}

FieldLayout
PointGrid::get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({COL,vector_tag,VL},{get_num_local_dofs(),vector_dim,nvl});
}

std::shared_ptr<AbstractGrid>
PointGrid::clone (const std::string& clone_name,
                  const bool shallow) const
{
  auto grid = std::make_shared<PointGrid> (clone_name,get_num_local_dofs(),get_num_vertical_levels(),get_comm());
  grid->copy_views(*this,shallow);
  return grid;
}

std::shared_ptr<PointGrid>
create_point_grid (const std::string& grid_name,
                   const int num_global_cols,
                   const int num_vertical_lev,
                   const ekat::Comm& comm)
{
  // Compute how many columns are owned by this rank
  const int num_procs = comm.size();

  auto num_my_cols = num_global_cols / num_procs;
  int remainder   = num_global_cols % num_procs;
  int dof_offset  = num_my_cols*comm.rank();
  if (comm.rank() < remainder) {
    ++num_my_cols;
    dof_offset += comm.rank();
  } else {
    dof_offset += remainder;
  }

  auto grid = std::make_shared<PointGrid>(grid_name,num_my_cols,num_vertical_lev,comm);
  grid->setSelfPointer(grid);

  PointGrid::dofs_list_type dofs_gids ("phys dofs",num_my_cols);
  auto h_dofs_gids = Kokkos::create_mirror_view(dofs_gids);
  std::iota(h_dofs_gids.data(),h_dofs_gids.data()+num_my_cols,dof_offset);
  Kokkos::deep_copy(dofs_gids,h_dofs_gids);

  grid->set_dofs(dofs_gids);

  using device_type       = DefaultDevice;
  using kokkos_types      = KokkosTypes<device_type>;
  using geo_view_type     = kokkos_types::view_1d<Real>;
  using C                 = scream::physics::Constants<Real>;

  // Store cell area, longitude, latitude and reference/surface pressure fractions
  // in geometry data. For  longitude, latitude and pressure fractions, set values
  // to NaN since they are not necessarily required from all application using PointGrid
  geo_view_type area("area", num_my_cols);
  geo_view_type lon ("lon",  num_my_cols);
  geo_view_type lat ("lat",  num_my_cols);
  geo_view_type hyam("hyam", num_vertical_lev);
  geo_view_type hybm("hybm", num_vertical_lev);

  // Estimate cell area for a uniform grid by taking the surface area
  // of the earth divided by the number of columns.  Note we do this in
  // units of radians-squared.
  const Real pi        = C::Pi;
  const Real cell_area = 4.0*pi/num_my_cols;

  Kokkos::deep_copy(area, cell_area);
  Kokkos::deep_copy(lon,  std::nan(""));
  Kokkos::deep_copy(lat,  std::nan(""));
  Kokkos::deep_copy(hyam, std::nan(""));
  Kokkos::deep_copy(hybm, std::nan(""));

  grid->set_geometry_data("area", area);
  grid->set_geometry_data("lon",  lon);
  grid->set_geometry_data("lat",  lat);
  grid->set_geometry_data("hyam", hyam);
  grid->set_geometry_data("hybm", hybm);

  return grid;
}

} // namespace scream
