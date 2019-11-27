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
  using string_pair = std::pair<std::string,std::string>;
  using remap_repo_type = std::map<string_pair,remapper_ptr_type>;

  UserProvidedGridsManager () = default;

  virtual ~UserProvidedGridsManager () = default;

  std::string name () const { return "User Provided Grids Manager"; }

  void set_remapper (const remapper_ptr_type remapper) {
    string_pair from_to = std::make_pair(remapper->get_src_grid()->name(),remapper->get_tgt_grid()->name());
    scream_require_msg (supports_grid(remapper->get_src_grid()->name()),
                        "Error! The remapper's source grid '" + from_to.first + "' is not supported."
                        "       Set the grids before setting the remappers.\n");
    scream_require_msg (supports_grid(remapper->get_tgt_grid()->name()),
                        "Error! The remapper's target grid '" + from_to.first + "' is not supported."
                        "       Set the grids before setting the remappers.\n");
    scream_require_msg (m_provided_remappers.find(from_to)==m_provided_remappers.end(),
                        "Error! A remapper from grid '" + from_to.first + "' to grid '" +
                        from_to.second + "' was already set.\n");
    m_provided_remappers[from_to] = remapper;
  }

  void set_grid (const grid_ptr_type grid) {
    scream_require_msg (m_provided_grids.find(grid->name())==m_provided_grids.end(),
                        "Error! A grid with name '" + grid->name() + "' was already set.\n");
    m_provided_grids[grid->name()] = grid;
  }

  void set_reference_grid (const std::string& grid_name) {
    scream_require_msg (m_provided_grids.find("Reference")==m_provided_grids.end(),
                        "Error! A reference was already set.\n");
    scream_require_msg (supports_grid(grid_name),
                        "Error! A grid with name '" + grid_name + "' was not set.\n");
    m_provided_grids["Reference"] = get_grid(grid_name);
  }

  void clean_up () {
    m_provided_remappers.clear();
    m_provided_grids.clear();
  }

protected:

  remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const {
    string_pair from_to = std::make_pair(from_grid->name(),to_grid->name());

    if (m_provided_remappers.find(from_to)!=m_provided_remappers.end()) {
      return m_provided_remappers.at(from_to);
    } else {
      scream_require_msg (false,
                          "Error! A remapper from grid '" + from_to.first + "' to grid '" +
                          from_to.second + "' (or viceversa) was not provided.\n");
    }
  }


  void build_grid (const std::string& grid_name) {
    // Simply make sure that the grid has been set
    scream_require_msg (supports_grid(grid_name),
                        "Error! No grid provided for '" + grid_name + "'.\n");
  }

  const grid_repo_type& get_repo () const { return m_provided_grids; }
        grid_repo_type& get_repo ()       { return m_provided_grids; }

  static grid_repo_type  m_provided_grids;

  static remap_repo_type m_provided_remappers;
};

inline std::shared_ptr<GridsManager>
create_user_provided_grids_manager (const Comm& /* comm */, const ParameterList& /* p */) {
  return std::make_shared<UserProvidedGridsManager>();
}

GridsManager::grid_repo_type UserProvidedGridsManager::m_provided_grids;
UserProvidedGridsManager::remap_repo_type UserProvidedGridsManager::m_provided_remappers;

} // namespace scream

#endif // SCREAM_USER_PROVIDED_GRIDS_MANAGER_HPP

