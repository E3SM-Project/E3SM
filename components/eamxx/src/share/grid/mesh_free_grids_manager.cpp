#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"
#include "share/property_checks/field_nan_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/io/scorpio_input.hpp"

#include "physics/share/physics_constants.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

#include <memory>
#include <numeric>

namespace scream {

MeshFreeGridsManager::
MeshFreeGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p)
 : m_params (p)
 , m_comm   (comm)
{
}

MeshFreeGridsManager::remapper_ptr_type
MeshFreeGridsManager::
do_create_remapper (const grid_ptr_type from_grid,
                    const grid_ptr_type to_grid) const
{
  return std::make_shared<DoNothingRemapper>(from_grid,to_grid);
}

void MeshFreeGridsManager::
build_grids ()
{
  auto has_positive_int = [&](const std::string& n) -> bool {
    return m_params.isParameter(n) && (m_params.get<int>(n)>0);
  };
  const bool build_pt = has_positive_int("number_of_global_columns");
  const bool build_se = has_positive_int("number_of_local_elements") &&
                        has_positive_int("number_of_gauss_points");

  const int num_vertical_levels = m_params.get<int>("number_of_vertical_levels");

  if (build_se) {
    // Build a set of completely disconnected spectral elements.
    const int num_local_elems  = m_params.get<int>("number_of_local_elements");
    const int num_gp           = m_params.get<int>("number_of_gauss_points");

    // Create the grid
    std::shared_ptr<SEGrid> se_grid;
    se_grid = std::make_shared<SEGrid>("SE Grid",num_local_elems,num_gp,num_vertical_levels,m_comm);
    se_grid->setSelfPointer(se_grid);

    // Set up the degrees of freedom.
    auto dof_gids = se_grid->get_dofs_gids();
    auto lid2idx  = se_grid->get_lid_to_idx_map();

    auto host_dofs    = dof_gids.template get_view<AbstractGrid::gid_type*,Host>();
    auto host_lid2idx = lid2idx.template get_view<int**,Host>();

    // Count unique local dofs. On all elems except the very last one (on rank N),
    // we have num_gp*(num_gp-1) unique dofs;
    int num_local_dofs = num_local_elems*num_gp*num_gp;
    int offset = num_local_dofs*m_comm.rank();

    for (int ie = 0; ie < num_local_elems; ++ie) {
      for (int igp = 0; igp < num_gp; ++igp) {
        for (int jgp = 0; jgp < num_gp; ++jgp) {
          int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
          int gid = offset + idof;
          host_dofs(idof) = gid;
          host_lid2idx(idof, 0) = ie;
          host_lid2idx(idof, 1) = igp;
          host_lid2idx(idof, 2) = jgp;
        }
      }
    }

    // Sync to device
    dof_gids.sync_to_dev();
    lid2idx.sync_to_dev();

    se_grid->m_short_name = "se";
    add_geo_data(se_grid);

    add_grid(se_grid);
  }
  if (build_pt) {
    const int num_global_cols  = m_params.get<int>("number_of_global_columns");
    auto pt_grid = create_point_grid("Point Grid",num_global_cols,num_vertical_levels,m_comm);

    const auto units = ekat::units::Units::nondimensional();

    auto area = pt_grid->create_geometry_data("area", pt_grid->get_2d_scalar_layout(), units);

    // Estimate cell area for a uniform grid by taking the surface area
    // of the earth divided by the number of columns.  Note we do this in
    // units of radians-squared.
    using PC             = scream::physics::Constants<Real>;
    const Real pi        = PC::Pi;
    const Real cell_area = 4.0*pi/num_global_cols;
    area.deep_copy(cell_area);
    area.sync_to_host();

    add_geo_data(pt_grid);
    pt_grid->m_short_name = "pt";

    add_grid(pt_grid);
    this->alias_grid("Point Grid", "Physics");
  }
}

void MeshFreeGridsManager::
add_geo_data (const nonconstgrid_ptr_type& grid) const
{
  if (!m_params.isParameter("geo_data_source")) {
    return;
  }

  // Load geo data fields if geo data filename is present, and file contains them
  const auto& geo_data_source = m_params.get<std::string>("geo_data_source");
  if (geo_data_source=="CREATE_EMPTY_DATA") {
    using namespace ShortFieldTagsNames;
    FieldLayout layout_mid ({LEV},{grid->get_num_vertical_levels()});
    const auto units = ekat::units::Units::nondimensional();

    auto lat  = grid->create_geometry_data("lat" , grid->get_2d_scalar_layout(), units);
    auto lon  = grid->create_geometry_data("lon" , grid->get_2d_scalar_layout(), units);
    auto hyam  = grid->create_geometry_data("hyam" , layout_mid, units);
    auto hybm  = grid->create_geometry_data("hybm" , layout_mid, units);

    lat.deep_copy(ekat::ScalarTraits<Real>::invalid());
    lon.deep_copy(ekat::ScalarTraits<Real>::invalid());
    hyam.deep_copy(ekat::ScalarTraits<Real>::invalid());
    hybm.deep_copy(ekat::ScalarTraits<Real>::invalid());
    lat.sync_to_dev();
    lon.sync_to_dev();
    hyam.sync_to_dev();
    hybm.sync_to_dev();
  } else if (geo_data_source=="IC_FILE"){
    const auto& filename = m_params.get<std::string>("ic_filename");
    if (scorpio::has_variable_c2f(filename.c_str(),"lat") &&
        scorpio::has_variable_c2f(filename.c_str(),"lon")) {
      load_lat_lon(grid,filename);
    }

    if (scorpio::has_variable_c2f(filename.c_str(),"hyam") &&
        scorpio::has_variable_c2f(filename.c_str(),"hybm")) {
      load_vertical_coordinates(grid,filename);
    }
  }
}

void MeshFreeGridsManager::
load_lat_lon (const nonconstgrid_ptr_type& grid, const std::string& filename) const
{
  using geo_view_host = AtmosphereInput::view_1d_host;
  const auto units = ekat::units::Units::nondimensional();

  auto lat  = grid->create_geometry_data("lat" , grid->get_2d_scalar_layout(), units);
  auto lon  = grid->create_geometry_data("lon" , grid->get_2d_scalar_layout(), units);

  // Create host mirrors for reading in data
  std::map<std::string,geo_view_host> host_views = {
    { "lat", lat.get_view<Real*,Host>() },
    { "lon", lon.get_view<Real*,Host>() }
  };

  // Store view layouts
  std::map<std::string,FieldLayout> layouts = {
    { "lat", lat.get_header().get_identifier().get_layout() },
    { "lon", lon.get_header().get_identifier().get_layout() }
  };

  // Read lat/lon into host views
  ekat::ParameterList lat_lon_reader_pl;
  lat_lon_reader_pl.set("Filename",filename);
  lat_lon_reader_pl.set<std::vector<std::string>>("Field Names",{"lat","lon"});

  AtmosphereInput lat_lon_reader(m_comm, lat_lon_reader_pl);
  lat_lon_reader.init(grid, host_views, layouts);
  lat_lon_reader.read_variables();
  lat_lon_reader.finalize();

  // Sync to dev
  lat.sync_to_dev();
  lon.sync_to_dev();

#ifndef NDEBUG
  for (auto f : {lat, lon}) {
  auto lat_check = std::make_shared<FieldNaNCheck>(lat,grid)->check();
  EKAT_REQUIRE_MSG (lat_check.result==CheckResult::Pass,
      "ERROR! NaN values detected in latitude field.\n" + lat_check.msg);
  auto lon_check = std::make_shared<FieldNaNCheck>(lon,grid)->check();
  EKAT_REQUIRE_MSG (lon_check.result==CheckResult::Pass,
      "ERROR! NaN values detected in longitude field.\n" + lon_check.msg);
  }
#endif
}

void MeshFreeGridsManager::
load_vertical_coordinates (const nonconstgrid_ptr_type& grid, const std::string& filename) const
{
  using geo_view_host = AtmosphereInput::view_1d_host;

  using namespace ShortFieldTagsNames;
  FieldLayout layout_mid ({LEV},{grid->get_num_vertical_levels()});
  const auto units = ekat::units::Units::nondimensional();

  auto hyam = grid->create_geometry_data("hyam", layout_mid, units);
  auto hybm = grid->create_geometry_data("hybm", layout_mid, units);

  // Create host mirrors for reading in data
  std::map<std::string,geo_view_host> host_views = {
    { "hyam", hyam.get_view<Real*,Host>() },
    { "hybm", hybm.get_view<Real*,Host>() }
  };

  // Store view layouts
  using namespace ShortFieldTagsNames;
  std::map<std::string,FieldLayout> layouts = {
    { "hyam", hyam.get_header().get_identifier().get_layout() },
    { "hybm", hybm.get_header().get_identifier().get_layout() }
  };

  // Read hyam/hybm into host views
  ekat::ParameterList vcoord_reader_pl;
  vcoord_reader_pl.set("Filename",filename);
  vcoord_reader_pl.set<std::vector<std::string>>("Field Names",{"hyam","hybm"});

  AtmosphereInput vcoord_reader(m_comm,vcoord_reader_pl);
  vcoord_reader.init(grid, host_views, layouts);
  vcoord_reader.read_variables();
  vcoord_reader.finalize();

  // Sync to dev
  hyam.sync_to_dev();
  hybm.sync_to_dev();

#ifndef NDEBUG
  for (auto f : {hyam, hybm}) {
    auto nan_check = std::make_shared<FieldNaNCheck>(f,grid)->check();
    EKAT_REQUIRE_MSG (nan_check.result==CheckResult::Pass,
        "ERROR! NaN values detected in " + f.name() + " field.\n" + nan_check.msg);
    auto interval_check = std::make_shared<FieldWithinIntervalCheck>(f,grid,0,1)->check();
    EKAT_REQUIRE_MSG (nan_check.result==CheckResult::Pass,
        "ERROR! Field " + f.name() + " has values not in [0,1].\n" + nan_check.msg);
  }
#endif
}

std::shared_ptr<GridsManager>
create_mesh_free_grids_manager (const ekat::Comm& comm, const int num_local_elems,
                                const int num_gp, const int num_vertical_levels,
                                const int num_global_cols)
{
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",num_global_cols);
  gm_params.set("number_of_local_elements",num_local_elems);
  gm_params.set("number_of_gauss_points",num_gp);
  gm_params.set("number_of_vertical_levels",num_vertical_levels);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  return gm;
}

} // namespace scream
