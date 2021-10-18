#include "share/grid//user_provided_grids_manager.hpp"

namespace scream {

GridsManager::grid_repo_type UserProvidedGridsManager::m_provided_grids;
UserProvidedGridsManager::remap_repo_type UserProvidedGridsManager::m_provided_remappers;
std::string UserProvidedGridsManager::m_ref_grid_name;

} // namespace scream
