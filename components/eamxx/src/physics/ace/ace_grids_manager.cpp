#include "ace_grids_manager.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "share/eamxx_types.hpp"
#include "share/grid/latlon_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"

namespace scream {

AceGridsManager::AceGridsManager(const ekat::Comm &comm,
                                 const ekat::ParameterList &p)
    : m_comm(comm), m_params(p) {
  // Nothing to do here
}

void AceGridsManager::build_grids() {
  auto nlat        = m_params.get<int>("num_latitude_points");
  auto nlon        = m_params.get<int>("num_longitude_points");
  auto nlev        = m_params.get<int>("num_vertical_levels");
  auto gaussian    = m_params.get<bool>("gaussian_grid");
  auto sc_ncolumns = m_params.get<int>("sc_ncolumns");
  auto sc_nlevels  = m_params.get<int>("sc_nlevels");

  EKAT_REQUIRE_MSG(
      gaussian,
      "Error! Currently ACE grids manager only supports gaussian_grid=true.\n");
  EKAT_REQUIRE_MSG(nlat > 0,
                   "Error! ACE grids manager requires positive value for "
                   "'num_latitude_points'.\n");
  EKAT_REQUIRE_MSG(nlon > 0,
                   "Error! ACE grids manager requires positive value for "
                   "'num_longitude_points'.\n");
  EKAT_REQUIRE_MSG(nlev > 0,
                   "Error! ACE grids manager requires positive value for "
                   "'num_vertical_levels'.\n");
  EKAT_REQUIRE_MSG(sc_ncolumns > 0,
                   "Error! ACE grids manager requires positive value for "
                   "'sc_ncolumns'.\n");
  EKAT_REQUIRE_MSG(sc_nlevels > 0,
                   "Error! ACE grids manager requires positive value for "
                   "'sc_nlevels'.\n");

  // Create *the* latlon grid
  auto latlon_grid = create_latlon_grid("AceLL", nlat, nlon, nlev, m_comm);
  add_nonconst_grid(latlon_grid);

  const auto nondim = ekat::units::Units::nondimensional();
  const ekat::units::Units deg(nondim, "deg");

  using namespace ShortFieldTagsNames;

  auto lat =
      latlon_grid->create_geometry_data("lat", FieldLayout({LAT}, {nlat}), deg);
  auto lon =
      latlon_grid->create_geometry_data("lon", FieldLayout({LON}, {nlon}), deg);
  auto area = latlon_grid->create_geometry_data(
      "area", FieldLayout({LAT, LON}, {nlat, nlon}), nondim);
  AtmosphereInput reader(m_params.get<std::string>("data_filename"),
                         latlon_grid, {lat, lon, area}, true);
  reader.read_variables();

  // Create PointGrid version of lat-lon grid
  auto pt_grid = create_point_grid("AcePG", nlat * nlon, nlev, m_comm);
  add_nonconst_grid(pt_grid);

  // Get lat/lon/area from lat-lon grid, and set it in the pt grid
  auto lat_pt = pt_grid->create_geometry_data(
      "lat", FieldLayout({COL}, {nlat * nlon}), deg);
  auto lon_pt = pt_grid->create_geometry_data(
      "lon", FieldLayout({COL}, {nlat * nlon}), deg);
  auto area_pt = pt_grid->create_geometry_data(
      "area", FieldLayout({COL}, {nlat * nlon}), nondim);
  auto lat_pt_h     = lat_pt.get_view<Real *, Host>();
  auto lon_pt_h     = lon_pt.get_view<Real *, Host>();
  auto area_pt_h    = area_pt.get_view<Real *, Host>();
  auto lat_ll_h     = lat.get_view<const Real *, Host>();
  auto lon_ll_h     = lon.get_view<const Real *, Host>();
  auto area_ll_h    = area.get_view<const Real **, Host>();
  auto lid_to_idx   = latlon_grid->get_lid_to_idx_map();
  auto lid_to_idx_h = lid_to_idx.get_view<const int **, Host>();
  for(int i = 0; i < pt_grid->get_num_local_dofs(); ++i) {
    auto ilat    = lid_to_idx_h(i, 0);
    auto ilon    = lid_to_idx_h(i, 1);
    lat_pt_h(i)  = lat_ll_h(ilat);
    lon_pt_h(i)  = lon_ll_h(ilon);
    area_pt_h(i) = area_ll_h(ilat, ilon);
  }
  lat_pt.sync_to_dev();
  lon_pt.sync_to_dev();
  area_pt.sync_to_dev();

  // Create a pg2 grid, assuming we always have to run on pg2
  auto sc_grid =
      create_point_grid("Physics PG2", sc_ncolumns, sc_nlevels, m_comm);

  // The cpl expects 1-based numbering for col gids
  auto sc_gids   = sc_grid->get_dofs_gids();
  auto sc_gids_h = sc_gids.get_view<int*,Host>();
  for (int icol=0; icol<sc_ncolumns; ++icol) {
    ++sc_gids_h(icol);
  }
  sc_gids.sync_to_dev();
  add_nonconst_grid(sc_grid);

  // Add an alias for the grid
  sc_grid->add_alias("Physics");

  auto lat_sc = sc_grid->create_geometry_data(
      "lat", FieldLayout({COL}, {sc_ncolumns}), deg);
  auto lon_sc = sc_grid->create_geometry_data(
      "lon", FieldLayout({COL}, {sc_ncolumns}), deg);
  auto area_sc = sc_grid->create_geometry_data(
      "area", FieldLayout({COL}, {sc_ncolumns}), nondim);
  auto frac_sc = sc_grid->create_geometry_data(
      "frac", FieldLayout({COL}, {sc_ncolumns}), nondim);
  auto mask_sc = sc_grid->create_geometry_data(
      "mask", FieldLayout({COL}, {sc_ncolumns}), nondim);
  AtmosphereInput reader_sc(m_params.get<std::string>("sc_data_filename"),
                            sc_grid, {lat_sc, lon_sc, area_sc, frac_sc, mask_sc}, true);
  reader_sc.read_variables();
}

}  // namespace scream
