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

  std::string name () const { return "User Provided Grids Manager"; }

  static void set_grid (const std::shared_ptr<grid_type> grid) {
    scream_require_msg (m_provided_grids.find(grid->name())==m_provided_grids.end(),
                        "Error! A grid with name '" + grid->name() + "' was already set.\n");
    m_provided_grids[grid->name()] = grid;
  }
protected:

  void build_grid (const std::string& grid_name) {
    // Simply make sure that the grid has been set
    scream_require_msg (supports_grid(grid_name),
                        "Error! No grid provided for '" + grid_name + "'.\n");
  }

  const repo_type& get_repo () const { return m_provided_grids; }

  // The static variable lets you put grids pointers in the manager even before
  // we dynamically create one.
  static repo_type  m_provided_grids;
};

inline std::shared_ptr<GridsManager>
create_user_provided_grids_manager (const Comm& /* comm */, const ParameterList& /* p */) {
  return std::make_shared<UserProvidedGridsManager>();
}

GridsManager::repo_type UserProvidedGridsManager::m_provided_grids;

} // namespace scream

#endif // SCREAM_USER_PROVIDED_GRIDS_MANAGER_HPP

