#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/scream_homme_interface.hpp"
#include "dynamics/homme/physics_dynamics_remapper.hpp"

#include "share/grid/se_grid.hpp"
#include "share/remap/inverse_remapper.hpp"

// Get all Homme's config properties
#include "hommexx_dimensions.hpp"

namespace scream
{

DynamicsDrivenGridsManager::DynamicsDrivenGridsManager (const Comm& /* comm */, const ParameterList& p)
 : m_params(p)
{
  // This is debatable: in theory, HommeDynamics should be crated *before*
  // this grids manager, so homme should already be inited.
  scream_require_msg (was_init_homme2_called_f90(),
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

  if (grid_name==m_params.get<std::string>("Reference Grid","Physics")) {
    m_grids["Reference"] = get_grid(grid_name);
  }
}

void DynamicsDrivenGridsManager::build_dynamics_grid () {
  if (m_grids.find("SE Dynamics")==m_grids.end()) {
    // Create dynamics dofs map
    const int nelemd = get_homme_param_value<int>("nelemd");
    SEGrid::dofs_map_type dyn_dofs("dyn dofs",nelemd*NP*NP);
    auto h_dyn_dofs = Kokkos::create_mirror_view(dyn_dofs);
    for (int ie=0; ie<nelemd; ++ie) {
      int gids[NP][NP];
      get_elem_cols_gids_f90 (ie, &gids[0][0]);

      for (int i=0; i<NP; ++i) {
        for (int j=0; j<NP; ++j) {
          auto dof_props = Kokkos::subview(h_dyn_dofs,ie*NP*NP+i*NP+j,Kokkos::ALL());
          dof_props(0) = ie;
          dof_props(1) = i;
          dof_props(2) = j;
          dof_props(3) = gids[i][j];
        }
      }
    }
    Kokkos::deep_copy(dyn_dofs,h_dyn_dofs);

    // Not much to do: simply create a DefaultGrid
    m_grids["SE Dynamics"] = std::make_shared<SEGrid>(dyn_dofs,"SE Dynamics",GridType::SE_CellBased);
  }
}

void DynamicsDrivenGridsManager::build_physics_grid () {
  if (m_grids.find("SE Physics")==m_grids.end()) {
    // Make sure the dynamics grid is built before the physics one.
    // If already built, this returns immediately.
    build_dynamics_grid();

    // Create the physics dofs map
    const int num_cols = get_num_owned_columns_f90 ();
    SEGrid::dofs_map_type phys_dofs("phys dofs",num_cols);
    auto h_phys_dofs = Kokkos::create_mirror_view(phys_dofs);
    std::set<int> gdofs;
    auto dyn_grid = std::dynamic_pointer_cast<SEGrid>(m_grids.at("SE Dynamics"));
    scream_require_msg(static_cast<bool>(dyn_grid),"Error! The dynamics grid is not of type SEGrid.\n");

    auto dyn_dofs = dyn_grid->get_dofs_map();
    auto h_dyn_dofs = Kokkos::create_mirror_view(dyn_dofs);
    Kokkos::deep_copy(h_dyn_dofs,dyn_dofs);
    for (int idof=0, k=0; idof<dyn_grid->get_num_dofs(); ++idof) {
      auto it_bool = gdofs.insert(h_dyn_dofs(idof,3));
      if (it_bool.second) {
        h_phys_dofs(k,0) = h_dyn_dofs(idof,0);
        h_phys_dofs(k,1) = h_dyn_dofs(idof,1);
        h_phys_dofs(k,2) = h_dyn_dofs(idof,2);
        h_phys_dofs(k,3) = h_dyn_dofs(idof,3);
        ++k;
      }
    }
    Kokkos::deep_copy(phys_dofs,h_phys_dofs);

    // Not much to do: simply create a DefaultGrid
    m_grids["SE Physics"] = std::make_shared<SEGrid>(phys_dofs,"SE Physics",GridType::SE_NodeBased);
  }
}

} // namespace scream
