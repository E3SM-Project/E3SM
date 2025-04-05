#include "share/grid/latlon_grid.hpp"

#include <numeric>

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {

LatLonGrid::LatLonGrid(const std::string &grid_name, const int nglat,
                       const int nglon, const int my_nlon,
                       const int num_vertical_levels, const ekat::Comm &comm)
    : AbstractGrid(grid_name, GridType::Point, nglat * my_nlon, nglat * nglon,
                   num_vertical_levels, comm) {
  m_nlat_global = nglat;
  m_nlon_local  = my_nlon;
  m_nlon_global = nglon;

  create_dof_fields(get_2d_scalar_layout().rank());

  int lon_offset = m_nlon_local;
  comm.scan(&lon_offset, 1, MPI_SUM);
  lon_offset -= m_nlon_local;  // scan is inclusive, but we need exclusive sum

  // The lid->idx map maps idof->(ilat,ilon)
  auto lid2idx      = get_lid_to_idx_map();
  auto h_lid_to_idx = lid2idx.get_view<int **, Host>();

  for(int ilat = 0, idof = 0; ilat < m_nlat_global; ++ilat) {
    for(int ilon = 0; ilon < m_nlon_local; ++ilon, ++idof) {
      h_lid_to_idx(idof, 0) = ilat;
      h_lid_to_idx(idof, 1) = ilon;
    }
  }
  lid2idx.sync_to_dev();
  using namespace ShortFieldTagsNames;

  // The partitioned dim is the longitude
  const auto units = ekat::units::Units::nondimensional();
  m_partitioned_dim_gids =
      Field(FieldIdentifier("lon_gids", FieldLayout({LON}, {m_nlon_local}),
                            units, this->name(), DataType::IntType));
  m_partitioned_dim_gids.allocate_view();
}

FieldLayout LatLonGrid::get_2d_scalar_layout() const {
  using namespace ShortFieldTagsNames;

  return FieldLayout({LAT, LON}, {m_nlat_global, m_nlon_local})
      .rename_dims(m_special_tag_names);
}

FieldLayout LatLonGrid::get_2d_vector_layout(
    const int vector_dim, const std::string &vec_dim_name) const {
  using namespace ShortFieldTagsNames;

  auto fl = get_2d_scalar_layout();
  fl.append_dim(CMP, vector_dim, vec_dim_name);
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout LatLonGrid::get_2d_tensor_layout(
    const std::vector<int> &cmp_dims,
    const std::vector<std::string> &cmp_names) const {
  EKAT_REQUIRE_MSG(
      cmp_names.size() == cmp_dims.size(),
      "[LatLonGrid::get_2d_tensor_layout] Input vector dimensions mismatch.\n"
      "  - grid name: " + name() + "\n"
      "  - cmp_names: " + ekat::join(cmp_names, ",") + "\n"
      "  - cmp_dims : " + ekat::join(cmp_dims, ",") + "\n");
  using namespace ShortFieldTagsNames;

  auto fl = get_2d_scalar_layout();

  for(size_t i = 0; i < cmp_dims.size(); ++i) {
    fl.append_dim(CMP, cmp_dims[i], cmp_names[i]);
  }
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout LatLonGrid::get_3d_scalar_layout(const bool midpoints) const {
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  auto fl = get_2d_scalar_layout();
  fl.append_dim(VL, nvl);
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout LatLonGrid::get_3d_vector_layout(
    const bool midpoints, const int vector_dim,
    const std::string &vec_dim_name) const {
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  auto fl = get_2d_scalar_layout();
  fl.append_dim(CMP, vector_dim, vec_dim_name);
  fl.append_dim(VL, nvl);
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout LatLonGrid::get_3d_tensor_layout(
    const bool midpoints, const std::vector<int> &cmp_dims,
    const std::vector<std::string> &cmp_names) const {
  EKAT_REQUIRE_MSG(
      cmp_names.size() == cmp_dims.size(),
      "[LatLonGrid::get_2d_tensor_layout] Input vector dimensions mismatch.\n"
      "  - grid name: " + name() + "\n"
      "  - cmp_names: " + ekat::join(cmp_names, ",") + "\n"
      "  - cmp_dims : " + ekat::join(cmp_dims, ",") + "\n");
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  auto fl = get_2d_scalar_layout();
  for(size_t i = 0; i < cmp_dims.size(); ++i) {
    fl.append_dim(CMP, cmp_dims[i], cmp_names[i]);
  }
  fl.append_dim(VL, nvl);

  return fl.rename_dims(m_special_tag_names);
}

std::shared_ptr<AbstractGrid> LatLonGrid::clone(const std::string &clone_name,
                                                const bool shallow) const {
  auto grid = std::make_shared<LatLonGrid>(
      clone_name, m_nlat_global, m_nlon_global, m_nlon_local,
      get_num_vertical_levels(), get_comm());
  grid->copy_data(*this, shallow);
  return grid;
}

std::shared_ptr<LatLonGrid> create_latlon_grid(const std::string &grid_name,
                                               const int nglat, const int nglon,
                                               const int num_vertical_lev,
                                               const ekat::Comm &comm) {
  // Divide evenly. If remainder R is present, add 1 for first R ranks
  int my_nlon = nglon / comm.size();
  int rem     = nglon % comm.size();
  if(comm.rank() < rem) ++my_nlon;

  auto grid = std::make_shared<LatLonGrid>(grid_name, nglat, nglon, my_nlon,
                                           num_vertical_lev, comm);
  grid->setSelfPointer(grid);

  return grid;
}

}  // namespace scream
