#include "physics/share/physics_only_grids_manager.hpp"

#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/se_grid.hpp"

namespace scream {

PhysicsOnlyGridsManager::
PhysicsOnlyGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p)
 : m_params (p)
 , m_comm   (comm)
{
  // Nothing else to do here
}

PhysicsOnlyGridsManager::remapper_ptr_type
PhysicsOnlyGridsManager::
do_create_remapper (const grid_ptr_type from_grid,
                    const grid_ptr_type to_grid) const
{
  // As of today (6/2020), we only support the old np4 grid for physis.
  // TODO: add support for pg2
  EKAT_REQUIRE_MSG(from_grid->name()==to_grid->name(),
                   "Error! So far, PhysicsOnlyGridsManager assumes only one type of grid for physiccs.\n");

  return std::make_shared<IdentityRemapper<Real,device_type>>(from_grid);
}

void PhysicsOnlyGridsManager::build_grid (const std::string& grid_name) {
  EKAT_REQUIRE_MSG (grid_name=="SE Physics" || grid_name=="Physics",
                    "Error! Only 'SE Physics' (and, for convenience, 'Physics') grid supported for physics-only runs.\n"
                    "       Requested grid: " + grid_name + "\n");

  const auto& phys_only_gm_params = m_params.sublist("Physics Only");
  const int num_global_cols = phys_only_gm_params.get<int>("Number of global columns");

  const int num_atm_procs = m_comm.size();

  int num_my_cols = num_global_cols / num_atm_procs;
  int remainder   = num_global_cols % num_atm_procs;
  int dof_offset  = num_my_cols*m_comm.rank();
  if (m_comm.rank() < remainder) {
    ++num_my_cols;
    dof_offset += m_comm.rank();
  } else {
    dof_offset += remainder;
  }

  SEGrid::dofs_list_type phys_dofs("phys dofs",num_my_cols);
  auto h_phys_dofs = Kokkos::create_mirror_view(phys_dofs);
  for (int i=0; i<num_my_cols; ++i) {
    h_phys_dofs(i) = dof_offset + i;
  }
  Kokkos::deep_copy(phys_dofs,h_phys_dofs);

  m_grids["SE Physics"] = std::make_shared<SEGrid>(phys_dofs,"SE Physics",GridType::SE_NodeBased);

  if (grid_name==m_params.get<std::string>("Reference Grid")) {
    m_grids["Reference"] = get_grid(grid_name);
  }

  // Allow the user to simply call this grid 'Physics'
  m_grids["Physics"] = m_grids.at("SE Physics");
}

} // namespace scream
