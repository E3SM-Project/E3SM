#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/scream_homme_interface.hpp"
#include "dynamics/homme/physics_dynamics_remapper.hpp"

#include "share/grid/se_grid.hpp"
#include "share/grid/remap/inverse_remapper.hpp"

// Get all Homme's config properties
#include "hommexx_dimensions.hpp"

namespace scream
{

DynamicsDrivenGridsManager::DynamicsDrivenGridsManager (const ekat::Comm& /* comm */, const ekat::ParameterList& p)
 : m_params(p)
{
  // This is debatable: in theory, HommeDynamics should be crated *before*
  // this grids manager, so homme should already be inited.
  EKAT_REQUIRE_MSG (was_init_homme2_called_f90(),
                      "Error! Homme needs to be inited before the DynamicsDrivenGridsManager is created.\n");
}

DynamicsDrivenGridsManager::remapper_ptr_type
DynamicsDrivenGridsManager::do_create_remapper (const grid_ptr_type from_grid,
                                                const grid_ptr_type to_grid) const {
  using PDR = PhysicsDynamicsRemapper<remapper_type::scalar_type,device_type>;
  auto pd_remapper = std::make_shared<PDR>(m_grids.at("SE Physics"),m_grids.at("SE Dynamics"));
  if (from_grid->name()=="SE Physics" &&
      to_grid->name()=="SE Dynamics") {
    return pd_remapper;
  } else if (from_grid->name()=="SE Dynamics" &&
             to_grid->name()=="SE Physics") {
    return std::make_shared<InverseRemapper<Real,device_type>>(pd_remapper);
  }
  return nullptr;
}

void DynamicsDrivenGridsManager::build_grid (const std::string& grid_name)
{
  if (grid_name=="SE Physics") {
    build_physics_grid();
  } else if (grid_name=="SE Dynamics") {
    build_dynamics_grid();
  }

  if (grid_name==m_params.get<std::string>("Reference Grid")) {
    m_grids["Reference"] = get_grid(grid_name);
  }
}

void DynamicsDrivenGridsManager::build_dynamics_grid () {
  if (m_grids.find("SE Dynamics")==m_grids.end()) {

    // Initialize the dyn grid
    const int nelemd = get_homme_param_value<int>("nelemd");
    auto dyn_grid = std::make_shared<SEGrid>("SE Dynamics",GridType::SE_CellBased,nelemd,NP);

    // Create dynamics dofs map
    SEGrid::dofs_list_type      dofs("dyn dofs",nelemd*NP*NP);
    SEGrid::lid_to_idx_map_type lids_to_elgpgp("dyn lid to elgpgp",nelemd*NP*NP,3);
    auto h_dofs = Kokkos::create_mirror_view(dofs);

    get_elem_cols_gids_f90(h_dofs.data());

    auto h_lids_to_elgpgp = Kokkos::create_mirror_view(lids_to_elgpgp);
    for (int ie=0; ie<nelemd; ++ie) {
      for (int i=0; i<NP; ++i) {
        for (int j=0; j<NP; ++j) {
          const int idof = ie*NP*NP+i*NP+j;
          h_lids_to_elgpgp(idof,0) = ie;
          h_lids_to_elgpgp(idof,1) = i;
          h_lids_to_elgpgp(idof,2) = j;
        }
      }
    }
    Kokkos::deep_copy(dofs,h_dofs);
    Kokkos::deep_copy(lids_to_elgpgp,h_lids_to_elgpgp);

    dyn_grid->set_dofs (dofs, lids_to_elgpgp);

    // Set the grid in the map
    m_grids["SE Dynamics"] = dyn_grid;
  }
}

void DynamicsDrivenGridsManager::build_physics_grid () {
  if (m_grids.find("SE Physics")==m_grids.end()) {

    // Initialize the phys grid
    const int nelemd = get_homme_param_value<int>("nelemd");
    auto phys_grid = std::make_shared<SEGrid>("SE Physics",GridType::SE_NodeBased,nelemd,NP);

    // Create the physics dofs map
    const int num_cols = get_num_owned_columns_f90 ();
    SEGrid::dofs_list_type      dofs("phys dofs",num_cols);
    SEGrid::lid_to_idx_map_type lids_to_elgpgp("phys lid to elgpgp", num_cols, 3);
    auto h_dofs = Kokkos::create_mirror_view(dofs);
    auto h_lids_to_elgpgp = Kokkos::create_mirror_view(lids_to_elgpgp);
    std::set<int> gdofs;

    // Get (ie,igp,jgp,gid) data for each dof
    get_unique_cols_f90(h_dofs.data(),h_lids_to_elgpgp.data());

    Kokkos::deep_copy(dofs,h_dofs);
    Kokkos::deep_copy(lids_to_elgpgp,h_lids_to_elgpgp);

    // Set the grid in the map
    m_grids["SE Physics"] = phys_grid;
  }
}

} // namespace scream
