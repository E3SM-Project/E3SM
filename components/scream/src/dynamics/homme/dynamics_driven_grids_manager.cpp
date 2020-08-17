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
    // Create dynamics dofs map
    const int nelemd = get_homme_param_value<int>("nelemd");
    SEGrid::dofs_list_type dyn_dofs("dyn dofs",nelemd*NP*NP);
    SEGrid::dofs_map_type dyn_dofs_map("dyn dofs",nelemd*NP*NP);
    auto h_dyn_dofs = Kokkos::create_mirror_view(dyn_dofs);

    get_elem_cols_gids_f90(h_dyn_dofs.data());

    auto h_dyn_dofs_map = Kokkos::create_mirror_view(dyn_dofs_map);
    for (int ie=0; ie<nelemd; ++ie) {
      for (int i=0; i<NP; ++i) {
        for (int j=0; j<NP; ++j) {
          const int idof = ie*NP*NP+i*NP+j;
          h_dyn_dofs_map(idof,0) = ie;
          h_dyn_dofs_map(idof,1) = i;
          h_dyn_dofs_map(idof,2) = j;
        }
      }
    }
    Kokkos::deep_copy(dyn_dofs,h_dyn_dofs);
    Kokkos::deep_copy(dyn_dofs_map,h_dyn_dofs_map);

    // Not much to do: simply create a DefaultGrid
    m_grids["SE Dynamics"] = std::make_shared<SEGrid>(dyn_dofs_map,dyn_dofs,"SE Dynamics",GridType::SE_CellBased);
  }
}

void DynamicsDrivenGridsManager::build_physics_grid () {
  if (m_grids.find("SE Physics")==m_grids.end()) {
    // Create the physics dofs map
    const int num_cols = get_num_owned_columns_f90 ();
    SEGrid::dofs_list_type phys_dofs("phys dofs",num_cols);
    SEGrid::dofs_map_type phys_dofs_map("phys dofs_map",num_cols);
    auto h_phys_dofs = Kokkos::create_mirror_view(phys_dofs);
    auto h_phys_dofs_map = Kokkos::create_mirror_view(phys_dofs_map);
    std::set<int> gdofs;

    // Get (ie,igp,jgp,gid) data for each dof
    get_unique_cols_f90(h_phys_dofs.data(),h_phys_dofs_map.data());

    Kokkos::deep_copy(phys_dofs,h_phys_dofs);
    Kokkos::deep_copy(phys_dofs_map,h_phys_dofs_map);

    // Not much to do: simply create a DefaultGrid
    m_grids["SE Physics"] = std::make_shared<SEGrid>(phys_dofs_map,phys_dofs,"SE Physics",GridType::SE_NodeBased);
  }
}

} // namespace scream
