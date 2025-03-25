#include "ace_grids_manager.hpp"

#include "share/io/eamxx_scorpio_interface.hpp"

namespace scream
{

AceGridsManager::
AceGridsManager (const ekat::Comm& comm,
                 const ekat::ParameterList& p)
 : m_comm (comm)
 , m_params(p)
{
  // Nothing to do here
}

void AceGridsManager::
build_grids ()
{
  auto nlat = m_params.get<int>("num_latitude_points");
  auto nlon = m_params.get<int>("num_longitude_points");
  auto nlev = m_params.get<int>("num_vertical_levels");
  auto gaussian = m_params.get<bool>("gaussian_grid");

  EKAT_REQUIRE_MSG (gaussian,
      "Error! Currently ACE grids manager only supports gaussian_grid=true.\n");
  EKAT_REQUIRE_MSG (nlat>0,
      "Error! ACE grids manager requires positive value for 'num_latitude_points'.\n");
  EKAT_REQUIRE_MSG (nlon>0,
      "Error! ACE grids manager requires positive value for 'num_longitude_points'.\n");
  EKAT_REQUIRE_MSG (nlev>0,
      "Error! ACE grids manager requires positive value for 'num_longitude_points'.\n");

  int my_nlon = nlon / m_comm.size();
  int rem = nlon % m_comm.size();
  if (m_comm.rank()<rem)
    ++my_nlon;

  // Create latlon grid
  auto latlon_grid = create_latlon_grid ("Physics LatLon",nlat,nlon,my_nlon,nlev,m_comm);
  add_nonconst_grid(latlon_grid);

  const auto nondim = ekat::units::Units::nondimensional();
  const ekat::units::Units deg(nondim,"deg");
  auto lat  = latlon_grid->create_geometry_data("lat" ,  latlon_grid->get_2d_scalar_layout(), deg);
  auto lon  = latlon_grid->create_geometry_data("lon" ,  latlon_grid->get_2d_scalar_layout(), deg);
  auto area = latlon_grid->create_geometry_data("area" , latlon_grid->get_2d_scalar_layout(), nondim);
  AtmosphereInput reader(m_params.get<std::string>("data_filename"),latlon_grid,{lat,lon,area},true);
  reader.read_variables();

  // Create PointGrid version of latlon grid
  auto pt_grid = create_latlon_grid("Physics",nlat*nlon,nlat*my_nlon,nlev,m_comm);
  add_nonconst_grid(pt_grid);

  // Get lat/lon/area from latlon grid, and set it in the pt grid
  auto lat_pt  = pt_grid->create_geometry_data("lat", pt_grid->get_2d_scalar_layout(),deg);
  auto lon_pt  = pt_grid->create_geometry_data("lon", pt_grid->get_2d_scalar_layout(),deg);
  auto area_pt = pt_grid->create_geometry_data("area",pt_grid->get_2d_scalar_layout(),nondim);
  auto lat_pt_h  = lat_pt.get_view<Real*,Host>();
  auto lon_pt_h  = lon_pt.get_view<Real*,Host>();
  auto area_pt_h = area_pt.get_view<Real*,Host>();
  auto lat_ll_h  = lat.get_view<const Real**,Host>();
  auto lon_ll_h  = lon.get_view<const Real**,Host>();
  auto area_ll_h = area.get_view<const Real**,Host>();
  auto lid_to_idx = latlon_grid->get_lid_to_idx_map();
  auto lid_to_idx_h = lid_to_idx.get_view<const int**,Host>();
  for (int i=0; i<pt_grid->get_num_local_dofs(); ++i) {
    auto ilat = lid_to_idx_h(i,0);
    auto ilon = lid_to_idx_h(i,1);

    lat_pt_h(i) = lat_ll_h(ilat,ilon);
    lon_pt_h(i) = lon_ll_h(ilat,ilon);
    area_pt_h(i) = area_ll_h(ilat,ilon);
  }
  lat_pt.sync_to_dev();
  lon_pt.sync_to_dev();
  area_pt.sync_to_dev();
}

} // namespace scream
