#include "dynamics/homme/homme_grids_manager.hpp"
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

HommeGridsManager::
HommeGridsManager (const ekat::Comm& comm,
                   const ekat::ParameterList& p)
 : m_comm (comm)
 , m_params(p)
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
    auto nlname = m_params.get<std::string>("dynamics_namelist_file_name").c_str();
    init_params_f90 (nlname);
  }

  // Check that the global number of 2d elements is no less than the number of MPI ranks
  EKAT_REQUIRE_MSG (get_homme_param<int>("nelem")>=comm.size(),
      "Error! We do not yet support running EAMxx with a number of MPI ranks\n"
      "       larger than the number of 2d elements in Homme.\n"
      "  - num MPI ranks: " + std::to_string(comm.size()) + "\n"
      "  - num 2d elems : " + std::to_string(get_homme_param<int>("nelem")) + "\n");

  // Create the grid integer codes map (i.e., int->string
  build_pg_codes ();
}

HommeGridsManager::
~HommeGridsManager () {
  // Cleanup the grids stuff
  finalize_geometry_f90 ();

  // This class is done with Homme. Remove from its users list
  HommeContextUser::singleton().remove_user();
}

HommeGridsManager::remapper_ptr_type
HommeGridsManager::do_create_remapper (const grid_ptr_type from_grid,
                                       const grid_ptr_type to_grid) const {
  const auto from = from_grid->name();
  const auto to   = to_grid->name();

  EKAT_REQUIRE_MSG ( from=="Dynamics" || to=="Dynamics",
    "Error! Either source or target grid must be 'Dynamics'.\n");

  const bool p2d = to=="Dynamics";

  if (from=="Physics GLL" || to=="Physics GLL") {
    using PDR = PhysicsDynamicsRemapper;

    auto dyn_grid = get_grid("Dynamics");
    auto phys_grid = get_grid("Physics GLL");

    auto pd_remapper = std::make_shared<PDR>(phys_grid,dyn_grid);
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

void HommeGridsManager::
build_grids ()
{
  // Get the physics grid specs
  const ci_string pg_type      = m_params.get<std::string>("physics_grid_type");
  const ci_string pg_rebalance = m_params.get<std::string>("physics_grid_rebalance","None");

  // Get the physics grid code
  std::vector<int> pg_codes {
    m_pg_codes["GLL"]["None"],  // We always need this to read/write dyn grid stuff
    m_pg_codes[pg_type][pg_rebalance]
  };
  // In case the two pg codes are the same...
  auto it = std::unique(pg_codes.begin(),pg_codes.end());
  const int* codes_ptr = pg_codes.data();
  init_grids_f90 (codes_ptr,std::distance(pg_codes.begin(),it));
  
  // We know we need the dyn grid, so build it
  build_dynamics_grid ();

  // Also the GLL grid with no rebalance is needed for sure
  build_physics_grid("GLL","None");

  // If (pg type,rebalance) is (GLL,None), this will be a no op
  build_physics_grid(pg_type,pg_rebalance);

  // Make "Physics" be an alias to whatever the pair (pg_type,rebalance) refers to
  std::string pg_name = "Physics " + pg_type;
  if (pg_rebalance!="None") {
    pg_name += " " + pg_rebalance;
  }

  this->alias_grid(pg_name,"Physics");

  // Clean up temporaries used during grid initialization
  cleanup_grid_init_data_f90 ();
}

void HommeGridsManager::build_dynamics_grid () {
  const std::string name = "Dynamics";
  if (has_grid(name)) {
    return;
  }

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

  add_grid(dyn_grid);
}

void HommeGridsManager::
build_physics_grid (const ci_string& type, const ci_string& rebalance) {
  std::string name = "Physics " + type;
  if (rebalance != "None") {
    name += " " + rebalance;
  }

  // Build only if not built yet
  if (has_grid(name)) {
    return;
  }

  // Get the grid pg_type
  const int pg_code = m_pg_codes.at(type).at(rebalance);

  // Get dimensions and create "empty" grid
  const int nlev  = get_nlev_f90();
  const int nlcols = get_num_local_columns_f90 (pg_code % 10);

  auto phys_grid = std::make_shared<PointGrid>(name,nlcols,nlev,m_comm);
  phys_grid->setSelfPointer(phys_grid);

  // Create the gids, coords, area views
  AbstractGrid::dofs_list_type dofs("phys dofs",nlcols);
  AbstractGrid::geo_view_type  lat("lat",nlcols);
  AbstractGrid::geo_view_type  lon("lon",nlcols);
  AbstractGrid::geo_view_type  area("area",nlcols);
  AbstractGrid::geo_view_type  hyam("hyam",nlev);
  AbstractGrid::geo_view_type  hybm("hybm",nlev);
  auto h_dofs = Kokkos::create_mirror_view(dofs);
  auto h_lat  = Kokkos::create_mirror_view(lat);
  auto h_lon  = Kokkos::create_mirror_view(lon);
  auto h_area = Kokkos::create_mirror_view(area);

  // For the following, set to NaN. They will need to be grabbed
  // from input files later.
  Kokkos::deep_copy(hyam, std::nan(""));
  Kokkos::deep_copy(hybm, std::nan(""));

  // Get all specs of phys grid cols (gids, coords, area)
  get_phys_grid_data_f90 (pg_code, h_dofs.data(), h_lat.data(), h_lon.data(), h_area.data());

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
  phys_grid->set_geometry_data("hyam",hyam);
  phys_grid->set_geometry_data("hybm",hybm);

  add_grid(phys_grid);
}

void HommeGridsManager::
build_pg_codes () {

  // Codes for the physics grids supported
  // Rules: a code has two digits, X and Y
  //  - X denotes the column rebalance choice:
  //    - 0: no rebalance
  //    - 1: twin columns rebalance
  //  - Y denotes the type of physics grid:
  //    - 0: GLL grid
  //    - N: FV phys grid, with NxN points

  m_pg_codes["GLL"]["None"] =  0;
  m_pg_codes["GLL"]["Twin"] = 10;
  m_pg_codes["PG2"]["None"] =  2;
  m_pg_codes["PG2"]["Twin"] = 12;
} 

} // namespace scream
