#include "remapper_manager.hpp"

namespace scream {

RemappersManager::
RemappersManager (const std::shared_ptr<const GridsManager>& gm)
 : m_grids_manager(gm)
{
  EKAT_REQUIRE_MSG (gm,
      "[RemappersManager] Error! Input grids manager pointer is invalid.\n");
}

auto RemappersManager::
create_remapper (const grid_ptr_type& /* from_grid */,
                const grid_ptr_type& /* to_grid */) const
 -> remapper_ptr_type
{
  EKAT_ERROR_MSG ("[RemappersManager] Error! Derived class did not override 'create_remapper'.\n");
}

auto RemappersManager::
create_remapper (const std::string& from_grid,
                 const std::string& to_grid) const
 -> remapper_ptr_type
{
  return create_remapper(m_grids_manager->get_grid(from_grid),
                         m_grids_manager->get_grid(to_grid));
}

} // namespace scream
