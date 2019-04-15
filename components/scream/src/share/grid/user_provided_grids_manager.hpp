#ifndef SCREAM_USER_PROVIDED_GRIDS_MANAGER_HPP
#define SCREAM_USER_PROVIDED_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

namespace scream
{

// This class is meant to be used for small unit tests, where we want to
// test the Atmosphere Driver (AD) capabilities, without bothering too much
// about grids-related features. This manager lets you set pre-built grids
// in it rather than building them inside the manager.
class UserProvidedGridsManager : public GridsManager
{
public:

  UserProvidedGridsManager () = default;

  virtual ~UserProvidedGridsManager () = default;

  void build_grids (const std::set<GridType>& grid_types) {
    // Simply make sure that all types have been set
    for (auto type : grid_types) {
      scream_require_msg (m_provided_grids.count(type)==1,
                          "Error! No grid provided for type '" + e2str(type) + "'.\n");
    }
  }

  static void set_grid (const std::shared_ptr<grid_type> grid) {
    m_provided_grids[grid->type()] = grid;
  }
protected:

  const repo_type& get_repo () const { return m_provided_grids; }

  // The static variable lets you put grids pointers in the manager even before
  // we dynamically create one.
  static repo_type  m_provided_grids;
};

inline GridsManager*
create_user_provided_grids_manager (const ParameterList& /* p */) {
  return new UserProvidedGridsManager();
}

GridsManager::repo_type UserProvidedGridsManager::m_provided_grids;

} // namespace scream

#endif // SCREAM_USER_PROVIDED_GRIDS_MANAGER_HPP

