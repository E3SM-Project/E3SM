#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/scream_homme_interface.hpp"

#include "share/grid/default_grid.hpp"

#ifdef HAVE_CONFIG_H
// Get all Homme's config properties
#include "config.h.c"
#endif

namespace scream
{

DynamicsDrivenGridsManager::DynamicsDrivenGridsManager (const Comm& /* comm */, const ParameterList& /* p */)
{
  // This is debatable: in theory, HommeDynamics should be crated *before*
  // this grids manager, so homme should already be inited.
  scream_require_msg (was_init_homme2_called_f90(), "Error! Homme needs to be inited before the DynamicsDrivenGridsManager is created.\n");

  // Create dynamics dofs map first
  const int nelemd = get_homme_param_value<int>("nelemd");
  m_dyn_dofs = decltype(m_dyn_dofs)("dyn dofs",nelemd*NP*NP);
  auto h_dyn_dofs  = Kokkos::create_mirror_view(m_dyn_dofs);
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
  Kokkos::deep_copy(m_dyn_dofs,h_dyn_dofs);
}

void DynamicsDrivenGridsManager::build_grid (const std::string& grid_name)
{
  if (grid_name=="Physics") {
    build_physics_grid();
  } else if (grid_name=="Dynamics") {
    build_dynamics_grid();
  }
}

void DynamicsDrivenGridsManager::build_dynamics_grid () {
  // Not much to do: simply create a DefaultGrid
  
  if (m_grids.find("Dynamics")==m_grids.end()) {
    m_grids["Dynamics"] = std::make_shared<DefaultGrid<GridType::Dynamics>>(m_dyn_dofs,"Dynamics");
  }
}

void DynamicsDrivenGridsManager::build_physics_grid () {
  // We have to get the number of owned columns, to properly size the phys dofs map
  int get_num_owned_columns_f90 ();

  if (m_grids.find("Dynamics")==m_grids.end()) {
    m_grids["Dynamics"] = std::make_shared<DefaultGrid<GridType::Dynamics>>(m_dyn_dofs,"Dynamics");
  }
}

} // namespace scream
