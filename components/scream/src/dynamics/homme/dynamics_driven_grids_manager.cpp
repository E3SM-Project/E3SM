#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"

#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/remap/inverse_remapper.hpp"

// Get all Homme's compile-time dims
#include "homme_dimensions.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream
{

DynamicsDrivenGridsManager::
DynamicsDrivenGridsManager (const ekat::Comm& comm,
                            const ekat::ParameterList& p)
 : m_comm (comm)
{
  if (!is_parallel_inited_f90()) {
    // While we're here, we can init homme's parallel session
    auto fcomm = MPI_Comm_c2f(comm.mpi_comm());
    init_parallel_f90 (fcomm);
  }

  // This class needs Homme's context, so register as a user
  HommeContextUser::singleton().add_user();

  if (!is_params_inited_f90()) {
    // While we're here, we can init homme's parameters
    auto nlname = p.get<std::string>("Dynamics Namelist File Name").c_str();
    init_params_f90 (nlname);
  }

  // Check that the global number of 2d elements is no less than the number of MPI ranks
  EKAT_REQUIRE_MSG (get_homme_param<int>("nelem")>=comm.size(),
      "Error! We do not yet support running EAMxx with a number of MPI ranks\n"
      "       larger than the number of 2d elements in Homme.\n"
      "  - num MPI ranks: " + std::to_string(comm.size()) + "\n"
      "  - num 2d elems : " + std::to_string(get_homme_param<int>("nelem")) + "\n");

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

  // Get the ref grid name
  m_ref_grid_name = p.get<std::string>("Reference Grid");
  EKAT_REQUIRE_MSG (ekat::contains(gn,m_ref_grid_name),
      "Error! Invalid reference grid name: " + m_ref_grid_name + "\n");

  // Create the grid integer codes map (i.e., int->string
  build_grid_codes ();
}

DynamicsDrivenGridsManager::
~DynamicsDrivenGridsManager () {
  // Cleanup the grids stuff
  finalize_geometry_f90 ();

  // This class is done with Homme. Remove from its users list
  HommeContextUser::singleton().remove_user();
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
    using PDR = PhysicsDynamicsRemapper;

    auto pd_remapper = std::make_shared<PDR>(m_grids.at("Physics GLL"),dyn_grid);
    if (p2d) {
      return pd_remapper;
    } else {
      return std::make_shared<InverseRemapper>(pd_remapper);
    }
  } else {
    ekat::error::runtime_abort("Error! P-D remapping only implemented for 'Physics GLL' phys grid.\n");
  }
  return nullptr;
}

void DynamicsDrivenGridsManager::
build_grids (const std::set<std::string>& grid_names)
{
  // Retrieve all physics grids codes
  std::vector<int> pg_codes;
  for (const auto& gn : grid_names) {
    // Sanity check first
    EKAT_REQUIRE_MSG (supported_grids().count(gn)==1,
                      "Error! Grid '" + gn + "' is not supported by this grid manager.\n");

    auto code = m_grid_codes.at(gn);
    // Dyn grid has a negative code, while phys grids codes are >=0.
    // We only need to store the phys grids codes.
    if (code>=0) {
      pg_codes.push_back(code);
    }
  }

  pg_codes.push_back(m_grid_codes.at(m_ref_grid_name));
  const int* codes_ptr = pg_codes.data();
  init_grids_f90 (codes_ptr,pg_codes.size());
  
  // We know we need the dyn grid, so build it
  build_dynamics_grid ();

  for (const auto& gn : grid_names) {
    if (gn!="Dynamics") {
      build_physics_grid(gn);
    }
  }

  if (m_grids.find(m_ref_grid_name)==m_grids.end()) {
    build_physics_grid(m_ref_grid_name);
  }

  if (m_grids.find("Physics GLL")==m_grids.end()) {
    build_physics_grid("Physics GLL");
  }

  // Clean up temporaries used during grid initialization
  cleanup_grid_init_data_f90 ();
}

void DynamicsDrivenGridsManager::build_dynamics_grid () {
  const std::string name = "Dynamics";
  if (m_grids.find(name)==m_grids.end()) {
    // Get dimensions and create "empty" grid
    const int nlelem = get_num_local_elems_f90();
    const int nlev   = get_nlev_f90();

    auto dyn_grid = std::make_shared<SEGrid>("Dynamics",nlelem,HOMMEXX_NP,nlev,m_comm);
    dyn_grid->setSelfPointer(dyn_grid);

    const int ndofs = nlelem*HOMMEXX_NP*HOMMEXX_NP;

    // Create the gids, elgpgp, coords, area views
    AbstractGrid::dofs_list_type      dg_dofs("dyn dofs",ndofs);
    AbstractGrid::dofs_list_type      cg_dofs("dyn dofs",ndofs);
    AbstractGrid::lid_to_idx_map_type elgpgp("dof idx",ndofs,3);
    AbstractGrid::geo_view_type       lat("lat",ndofs);
    AbstractGrid::geo_view_type       lon("lon",ndofs);
    auto h_cg_dofs   = Kokkos::create_mirror_view(cg_dofs);
    auto h_dg_dofs   = Kokkos::create_mirror_view(dg_dofs);
    auto h_elgpgp = Kokkos::create_mirror_view(elgpgp);
    auto h_lat    = Kokkos::create_mirror_view(lat);
    auto h_lon    = Kokkos::create_mirror_view(lon);

    // Get (ie,igp,jgp,gid) data for each dof
    get_dyn_grid_data_f90 (h_dg_dofs.data(),h_cg_dofs.data(),h_elgpgp.data(), h_lat.data(), h_lon.data());

    Kokkos::deep_copy(dg_dofs,h_dg_dofs);
    Kokkos::deep_copy(cg_dofs,h_cg_dofs);
    Kokkos::deep_copy(elgpgp,h_elgpgp);
    Kokkos::deep_copy(lat,h_lat);
    Kokkos::deep_copy(lon,h_lon);

#ifndef NDEBUG
    // Check that latitude and longitude are valid
    using KT = KokkosTypes<DefaultDevice>;
    const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ndofs,nlev);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
      const Int i = team.league_rank();

      EKAT_KERNEL_ASSERT_MSG(!ekat::impl::is_nan(lat(i)), "Error! NaN values detected for latitude.");
      EKAT_KERNEL_ASSERT_MSG(!ekat::impl::is_nan(lon(i)), "Error! NaN values detected for longitude.");
    });
#endif

    // Set dofs and geo data in the grid
    dyn_grid->set_dofs (dg_dofs);
    dyn_grid->set_cg_dofs (cg_dofs);
    dyn_grid->set_lid_to_idx_map(elgpgp);
    dyn_grid->set_geometry_data ("lat", lat);
    dyn_grid->set_geometry_data ("lon", lon);

    m_grids["Dynamics"] = dyn_grid;
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
    const int nlcols = get_num_local_columns_f90 (pg_type % 10);

    auto phys_grid = std::make_shared<PointGrid>(name,nlcols,nlev,m_comm);
    phys_grid->setSelfPointer(phys_grid);

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
    
#ifndef NDEBUG
    // Check that latitude, longitude, and area are valid
    using KT = KokkosTypes<DefaultDevice>;
    const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(nlcols,nlev);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
      const Int i = team.league_rank();

      EKAT_KERNEL_ASSERT_MSG(!ekat::impl::is_nan(lat(i)),                   "Error! NaN values detected for latitude.");
      EKAT_KERNEL_ASSERT_MSG(!ekat::impl::is_nan(lon(i)),                   "Error! NaN values detected for longitude.");
      EKAT_KERNEL_ASSERT_MSG(area(i) >  0  && !ekat::impl::is_nan(area(i)), "Error! Non-positve or NaN values detected for area.");
    });
#endif

    // Set dofs and geo data in the grid
    phys_grid->set_dofs(dofs);
    phys_grid->set_geometry_data("lat",lat);
    phys_grid->set_geometry_data("lon",lon);
    phys_grid->set_geometry_data("area",area);

    m_grids[name] = phys_grid;
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
  constexpr int pg4_t = 14;  // Physics PG4 Twin

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
