#include "physics/share/physics_only_grids_manager.hpp"

#include "share/grid/simple_grid.hpp"
#include "share/grid/remap/identity_remapper.hpp"

namespace scream {
namespace physics {

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
  EKAT_REQUIRE_MSG (grid_name=="Physics Only" || grid_name=="Physics",
                    "Error! Only 'Physics Only' (and, for convenience, 'Physics') grid supported for physics-only runs.\n"
                    "       Requested grid: " + grid_name + "\n");

  const auto& phys_only_gm_params = m_params.sublist("Physics Only");
  const int num_global_cols = phys_only_gm_params.get<int>("Number of global columns");
  const int num_vertical_lev = phys_only_gm_params.get<int>("Number of vertical levels");

  auto grid = std::make_shared<SimpleGrid>("Physics",num_global_cols,num_vertical_lev,m_comm);
  m_grids["Physics Only"] = m_grids["Physics"] = grid;

  if (grid_name==m_params.get<std::string>("Reference Grid")) {
    m_grids["Reference"] = get_grid(grid_name);
  }
}

} // namespace physics
} // namespace scream
