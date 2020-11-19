#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "dynamics/homme/physics_dynamics_remapper.hpp"

#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/remap/inverse_remapper.hpp"

// Get all Homme's compile-time dims
#include "homme_dimensions.hpp"

namespace scream
{

DynamicsDrivenGridsManager::
DynamicsDrivenGridsManager (const ekat::Comm& /* comm */,
                            const ekat::ParameterList& /* p */)
{
  // Valid names for the dyn grid
  auto& gn = m_valid_grid_names;

  // Dynamics
  gn.insert("Dynamics");

  // Physics on same element-partition as dynamics
  gn.insert("Physics GLL"); // Phys columns are SE gll points
  gn.insert("Physics PG2"); // Phys columns are FV points in 2x2 subcells of SE cell
  gn.insert("Physics PG3"); // Phys columns are FV points in 3x3 subcells of SE cell
  gn.insert("Physics PG4"); // Phys columns are FV points in 4x4 subcells of SE cell

  // Physics with columns partitioned so that each rank owns "twin" columns
  gn.insert("Physics GLL Twin"); // Phys columns are SE gll points
  gn.insert("Physics PG2 Twin"); // Phys columns are FV points in 2x2 subcells of SE cell
  gn.insert("Physics PG3 Twin"); // Phys columns are FV points in 3x3 subcells of SE cell
  gn.insert("Physics PG4 Twin"); // Phys columns are FV points in 4x4 subcells of SE cell

  // TODO: add other rebalancing?

  // Create the grid integer codes map (i.e., int->string
  build_grid_codes ();
}

DynamicsDrivenGridsManager::
~DynamicsDrivenGridsManager () {
  // Cleanup the grids stuff
  finalize_geometry_f90 ();
}

DynamicsDrivenGridsManager::remapper_ptr_type
DynamicsDrivenGridsManager::do_create_remapper (const grid_ptr_type from_grid,
                                                const grid_ptr_type to_grid) const {
  const auto from = from_grid->name();
  const auto to   = to_grid->name();

  EKAT_REQUIRE_MSG ( from=="Dynamics" || to=="Dynamics",
    "Error! Either source or target grid must be 'Dynamics'.\n");

  const bool p2d = to=="Dynamics";

  auto dyn_grid = m_grids.at("Dynamics");

  if (from=="Physics GLL" || to=="Physics GLL") {
    using PDR = PhysicsDynamicsRemapper<remapper_type::real_type>;

    auto pd_remapper = std::make_shared<PDR>(m_grids.at("Physics"),dyn_grid);
    if (p2d) {
      return pd_remapper;
    } else {
      return std::make_shared<InverseRemapper<Real>>(pd_remapper);
    }
  } else {
    ekat::error::runtime_abort("Error! P-D remapping only implemented for 'Physics GLL' phys grid.\n");
  }
  return nullptr;
}

void DynamicsDrivenGridsManager::
build_grids (const std::set<std::string>& grid_names,
             const std::string& reference_grid) {
  // Retrieve all grid codes
  std::set<int> codes;
  for (const auto& gn : grid_names) {
    // Sanity check first
    EKAT_REQUIRE_MSG (supported_grids().count(gn)==1,
                      "Error! Grid '" + gn + "' is not supported by this grid manager.\n");

    codes.insert(m_grid_codes.at(gn));
  }

  // Deduce the phys grid type we need
  int pgN = -1;
  for (auto code : codes) {
    if (code>=0) {
      int N = code % 10;
      EKAT_REQUIRE_MSG (pgN==-1 || N==pgN,
                        "Error! Mixing different types of physics grid is not allowed.\n"
                        "       You can, however, have phys grids with different *balancing* options.\n");
      pgN = N;
    }
  }

  // Nobody should have init-ed the geometries yet. So error out if someone did.
  EKAT_REQUIRE_MSG (!is_geometry_inited_f90(), "Error! Geometry was somehow already init-ed.\n");

  init_grids_f90 (pgN);

  // We know we need the dyn grid, so build it
  build_dynamics_grid ();

  for (const auto& gn : grid_names) {
    if (gn!="Dynamics") {
      build_physics_grid(gn);
    }
  }

  // Now we can cleanup all the grid stuff in f90
  cleanup_grid_init_data_f90 ();

  // Set the ptr to the ref grid
  m_grids["Reference"] = get_grid(reference_grid); 
}

void DynamicsDrivenGridsManager::build_dynamics_grid () {
  if (m_grids.find("Dynamics")==m_grids.end()) {

    // Get dimensions and create "empty" grid
    const int nlelem = get_num_local_elems_f90();
    const int ngelem = get_num_global_elems_f90();
    const int nlev   = get_nlev_f90();
    auto dyn_grid = std::make_shared<SEGrid>("Dynamics",ngelem,nlelem,NP,nlev);

    const int ndofs = nlelem*NP*NP;

    // Create the gids, elgpgp, coords, area views
    AbstractGrid::dofs_list_type      dofs("dyn dofs",ndofs);
    AbstractGrid::lid_to_idx_map_type elgpgp("dof idx",ndofs,3);
    AbstractGrid::geo_view_type       lat("lat",ndofs);
    AbstractGrid::geo_view_type       lon("lon",ndofs);
    auto h_dofs   = Kokkos::create_mirror_view(dofs);
    auto h_elgpgp = Kokkos::create_mirror_view(elgpgp);
    auto h_lat    = Kokkos::create_mirror_view(lat);
    auto h_lon    = Kokkos::create_mirror_view(lon);

    // Get (ie,igp,jgp,gid) data for each dof
    get_dyn_grid_data_f90 (h_dofs.data(),h_elgpgp.data(), h_lat.data(), h_lon.data());

    Kokkos::deep_copy(dofs,h_dofs);
    Kokkos::deep_copy(elgpgp,h_elgpgp);
    Kokkos::deep_copy(lat,h_lat);
    Kokkos::deep_copy(lon,h_lon);

    // Set dofs and geo data in the grid
    dyn_grid->set_dofs (dofs, elgpgp);
    dyn_grid->set_geometry_data ("lat", lat);
    dyn_grid->set_geometry_data ("lon", lon);
  }
}

void DynamicsDrivenGridsManager::
build_physics_grid (const std::string& name) {

  // Build only if not built yet
  if (m_grids.find(name)==m_grids.end()) {

    // Get the grid pg_type
    const int pg_type = m_grid_codes.at(name);

    // Get dimensions and create "empty" grid
    const int nlev  = get_nlev_f90();
    const int nlcols = get_num_local_columns_f90 ();
    const int ngcols = get_num_global_columns_f90 ();

    auto phys_grid = std::make_shared<PointGrid>("Physics",ngcols,nlcols,nlev);

    // Create the gids, coords, area views
    AbstractGrid::dofs_list_type dofs("phys dofs",nlcols);
    AbstractGrid::geo_view_type  lat("lat",nlcols);
    AbstractGrid::geo_view_type  lon("lon",nlcols);
    AbstractGrid::geo_view_type  area("area",nlcols);
    auto h_dofs = Kokkos::create_mirror_view(dofs);
    auto h_lat  = Kokkos::create_mirror_view(lat);
    auto h_lon  = Kokkos::create_mirror_view(lon);
    auto h_area = Kokkos::create_mirror_view(area);

    // Get all specs of phys grid cols (gids, coords, area)
    get_phys_grid_data_f90 (pg_type, h_dofs.data(), h_lat.data(), h_lon.data(), h_area.data());

    Kokkos::deep_copy(dofs,h_dofs);
    Kokkos::deep_copy(lat, h_lat);
    Kokkos::deep_copy(lon, h_lon);
    Kokkos::deep_copy(area,h_area);

    // Set dofs and geo data in the grid
    phys_grid->set_dofs(dofs);
    phys_grid->set_geometry_data("lat",lat);
    phys_grid->set_geometry_data("lon",lon);
    phys_grid->set_geometry_data("area",area);
  }
}

void DynamicsDrivenGridsManager::
build_grid_codes () {

  // Codes for the physics grids supported
  // Rules:
  //  -  -1 is dynamics
  //  -  XY is physics, with X=1 denoting twin columns partition,
  //     X=0 denoting no rebalance, and Y being the value of N
  //     in the FV phys grid (N=2,3,4), with N=0 denoting GLL grid.
  constexpr int dyn   = -1;  // Dyanamics grid (not a phys grid)
  constexpr int gll   =  0;  // Physics GLL
  constexpr int pg2   =  2;  // Physics PG2
  constexpr int pg3   =  3;  // Physics PG3
  constexpr int pg4   =  3;  // Physics PG4
  constexpr int gll_t = 10;  // Physics GLL Twin
  constexpr int pg2_t = 12;  // Physics PG2 Twin
  constexpr int pg3_t = 13;  // Physics PG3 Twin
  constexpr int pg4_t = 13;  // Physics PG4 Twin

  for (const auto& name : m_valid_grid_names) {
    int code;

    if (name=="Physics GLL") {
      code = gll;
    } else if (name=="Physics GLL Twin") {
      code = gll_t;
    } else if (name=="Physics PG2") {
      code = pg2;
    } else if (name=="Physics PG2 Twin") {
      code = pg2_t;
    } else if (name=="Physics PG3") {
      code = pg3;
    } else if (name=="Physics PG3 Twin") {
      code = pg3_t;
    } else if (name=="Physics PG4") {
      code = pg4;
    } else if (name=="Physics PG4 Twin") {
      code = pg4_t;
    } else {
      code = dyn;
    }

    m_grid_codes[name] = code;
  } 
} 

} // namespace scream
