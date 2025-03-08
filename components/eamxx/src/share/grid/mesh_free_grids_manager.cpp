#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/property_checks/field_nan_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
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
  // Nothing else to do here
}

MeshFreeGridsManager::remapper_ptr_type
MeshFreeGridsManager::
do_create_remapper (const grid_ptr_type /* from_grid */,
                    const grid_ptr_type /* to_grid */) const
{
  EKAT_ERROR_MSG ("Error! MeshFreeGridsManager does not offer any remapper.\n");
  return nullptr;
}

void MeshFreeGridsManager::
build_grids ()
{
  const auto& names = m_params.get<std::vector<std::string>>("grids_names");
  for (const auto& gname : names) {
    auto& params = m_params.sublist(gname);
    const auto& type = params.get<std::string>("type");
    if (type=="se_grid") {
      build_se_grid (gname,params);
    } else if (type=="point_grid") {
      build_point_grid (gname,params);
    } else {
      EKAT_ERROR_MSG (
          "[MeshFreeGridsManager::build_grids] Unrecognized grid type.\n"
          " - grid name: " + gname + "\n"
          " - grid type: " + type + "\n"
          " - valid types: se_grid, point_grid\n");
    }
    const auto& aliases = params.get<std::vector<std::string>>("aliases",{});
    for (const auto& n : aliases) {
      this->alias_grid(gname, n);
    }
  }
}

void MeshFreeGridsManager::
build_se_grid (const std::string& name, ekat::ParameterList& params)
{
  // Build a set of completely disconnected spectral elements.
  const int num_local_elems  = params.get<int>("number_of_local_elements");
  const int num_gp           = params.get<int>("number_of_gauss_points");
  const int num_vertical_levels = params.get<int>("number_of_vertical_levels");

  // Create the grid
  std::shared_ptr<SEGrid> se_grid;
  se_grid = std::make_shared<SEGrid>(name,num_local_elems,num_gp,num_vertical_levels,m_comm);
  se_grid->setSelfPointer(se_grid);

  // Set up the degrees of freedom.
  auto dof_gids  = se_grid->get_dofs_gids();
  auto elem_gids = se_grid->get_partitioned_dim_gids();
  auto lid2idx   = se_grid->get_lid_to_idx_map();

  auto host_dofs    = dof_gids.template get_view<AbstractGrid::gid_type*,Host>();
  auto host_elems   = elem_gids.template get_view<AbstractGrid::gid_type*,Host>();
  auto host_lid2idx = lid2idx.template get_view<int**,Host>();

  // Count unique local dofs. On all elems except the very last one (on rank N),
  // we have num_gp*(num_gp-1) unique dofs;
  int num_local_dofs = num_local_elems*num_gp*num_gp;
  int offset = num_local_dofs*m_comm.rank();

  for (int ie = 0; ie < num_local_elems; ++ie) {
    host_elems[ie] = ie + num_local_elems*m_comm.rank();
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
  elem_gids.sync_to_dev();
  lid2idx.sync_to_dev();

  se_grid->m_short_name = "se";
  add_geo_data(se_grid);

  add_nonconst_grid(se_grid);
}

void MeshFreeGridsManager::
build_point_grid (const std::string& name, ekat::ParameterList& params)
{
  const int num_global_cols = params.get<int>("number_of_global_columns");
  const int num_vertical_levels = params.get<int>("number_of_vertical_levels");
  auto pt_grid = create_point_grid(name,num_global_cols,num_vertical_levels,m_comm);

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

  add_nonconst_grid(pt_grid);
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
    FieldLayout layout_int ({ILEV},{grid->get_num_vertical_levels()+1});
    const auto units = ekat::units::Units::nondimensional();

    auto lat  = grid->create_geometry_data("lat" ,  grid->get_2d_scalar_layout(), units);
    auto lon  = grid->create_geometry_data("lon" ,  grid->get_2d_scalar_layout(), units);
    auto hyam = grid->create_geometry_data("hyam" , layout_mid, units);
    auto hybm = grid->create_geometry_data("hybm" , layout_mid, units);
    auto hyai = grid->create_geometry_data("hyai" , layout_int, units);
    auto hybi = grid->create_geometry_data("hybi" , layout_int, units);
    auto lev  = grid->create_geometry_data("lev" ,  layout_mid, units);
    auto ilev = grid->create_geometry_data("ilev" , layout_int, units);

    lat.deep_copy(ekat::ScalarTraits<Real>::invalid());
    lon.deep_copy(ekat::ScalarTraits<Real>::invalid());
    hyam.deep_copy(ekat::ScalarTraits<Real>::invalid());
    hybm.deep_copy(ekat::ScalarTraits<Real>::invalid());
    lev.deep_copy(ekat::ScalarTraits<Real>::invalid());
    ilev.deep_copy(ekat::ScalarTraits<Real>::invalid());
    lat.sync_to_dev();
    lon.sync_to_dev();
    hyam.sync_to_dev();
    hybm.sync_to_dev();
    lev.sync_to_dev();
    ilev.sync_to_dev();
  } else if (geo_data_source=="IC_FILE"){
    const auto& filename = m_params.get<std::string>("ic_filename");
    if (scorpio::has_var(filename,"lat") &&
        scorpio::has_var(filename,"lon")) {
      load_lat_lon(grid,filename);
    }

    if (scorpio::has_var(filename,"hyam") &&
        scorpio::has_var(filename,"hybm") &&
        scorpio::has_var(filename,"hyai") &&
        scorpio::has_var(filename,"hybi") ) {
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

  AtmosphereInput lat_lon_reader(lat_lon_reader_pl, grid, host_views, layouts);
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
  using namespace ekat::units;

  FieldLayout layout_mid ({LEV},{grid->get_num_vertical_levels()});
  FieldLayout layout_int ({ILEV},{grid->get_num_vertical_levels()+1});
  Units nondim = Units::nondimensional();
  Units mbar (100*Pa,"mb");

  auto hyam = grid->create_geometry_data("hyam", layout_mid, nondim);
  auto hybm = grid->create_geometry_data("hybm", layout_mid, nondim);
  auto hyai = grid->create_geometry_data("hyai", layout_int, nondim);
  auto hybi = grid->create_geometry_data("hybi", layout_int, nondim);
  auto lev  = grid->create_geometry_data("lev",  layout_mid, mbar);
  auto ilev = grid->create_geometry_data("ilev", layout_int, mbar);

  // Create host mirrors for reading in data
  std::map<std::string,geo_view_host> host_views = {
    { "hyam", hyam.get_view<Real*,Host>() },
    { "hybm", hybm.get_view<Real*,Host>() },
    { "hyai", hyai.get_view<Real*,Host>() },
    { "hybi", hybi.get_view<Real*,Host>() }
  };

  // Store view layouts
  using namespace ShortFieldTagsNames;
  std::map<std::string,FieldLayout> layouts = {
    { "hyam", hyam.get_header().get_identifier().get_layout() },
    { "hybm", hybm.get_header().get_identifier().get_layout() },
    { "hyai", hyai.get_header().get_identifier().get_layout() },
    { "hybi", hybi.get_header().get_identifier().get_layout() }
  };

  // Read hyam/hybm into host views
  ekat::ParameterList vcoord_reader_pl;
  vcoord_reader_pl.set("Filename",filename);
  vcoord_reader_pl.set<std::vector<std::string>>("Field Names",{"hyam","hybm","hyai","hybi"});

  AtmosphereInput vcoord_reader(vcoord_reader_pl,grid, host_views, layouts);
  vcoord_reader.read_variables();
  vcoord_reader.finalize();

  // Build lev and ilev from hyam and hybm, and ilev from hyai and hybi
  using PC             = scream::physics::Constants<Real>;
  const Real ps0        = PC::P0;

  auto hyam_v  = hyam.get_view<const Real*,Host>();
  auto hybm_v  = hybm.get_view<const Real*,Host>();
  auto lev_v   = lev.get_view<Real*,Host>();
  auto hyai_v  = hyai.get_view<const Real*,Host>();
  auto hybi_v  = hybi.get_view<const Real*,Host>();
  auto ilev_v  = ilev.get_view<Real*,Host>();
  auto num_lev = grid->get_num_vertical_levels();
  for (int ii=0;ii<num_lev;ii++) {
    lev_v(ii)  = 0.01*ps0*(hyam_v(ii)+hybm_v(ii));
    ilev_v(ii) = 0.01*ps0*(hyai_v(ii)+hybi_v(ii));
  }
  // Note, ilev is just 1 more level than the number of midpoint levs
  ilev_v(num_lev) = 0.01*ps0*(hyai_v(num_lev)+hybi_v(num_lev));

  // Sync to dev
  hyam.sync_to_dev();
  hybm.sync_to_dev();
  hyai.sync_to_dev();
  hybi.sync_to_dev();
  lev.sync_to_dev();
  ilev.sync_to_dev();

#ifndef NDEBUG
  for (auto f : {hyam, hybm, hyai, hybi}) {
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
  std::vector<std::string> grids_names;
  if (num_local_elems>=2) {
    grids_names.push_back("SE Grid");
    auto& pl = gm_params.sublist("SE Grid");
    pl.set("type",std::string("se_grid"));
    pl.set("number_of_local_elements",num_local_elems);
    pl.set("number_of_gauss_points",num_gp);
    pl.set("number_of_vertical_levels",num_vertical_levels);
  }
  if (num_global_cols>0) {
    grids_names.push_back("Point Grid");
    auto& pl = gm_params.sublist("Point Grid");
    pl.set("type",std::string("point_grid"));
    pl.set("number_of_global_columns",num_global_cols);
    pl.set("number_of_vertical_levels",num_vertical_levels);
  }
  gm_params.set("grids_names",grids_names);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  return gm;
}

} // namespace scream
