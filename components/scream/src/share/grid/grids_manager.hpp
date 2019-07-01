#ifndef SCREAM_GRIDS_MANAGER_HPP
#define SCREAM_GRIDS_MANAGER_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_assert.hpp"
#include "share/util/factory.hpp"
#include "share/util/string_utils.hpp"
#include <map>
#include <set>
#include <memory>

namespace scream
{

class GridsManager
{
public:
  using grid_type     = AbstractGrid;
  using grid_ptr_type = std::shared_ptr<grid_type>;
  using repo_type     = std::map<std::string, grid_ptr_type>;

  GridsManager () = default;
  virtual ~GridsManager () = default;

  virtual std::string name () const = 0;

  std::shared_ptr<grid_type> get_grid (const std::string& name) const;

  void build_grids (const std::set<std::string>& grid_names) {
    for (const auto& name : grid_names) {
      build_grid(name);
    }
  }

protected:

  bool supports_grid (const std::string& grid_name) const {
    const auto& grids = get_repo ();
    return grids.find(grid_name)!=grids.end();
  }

  // This mini-function simply prints the supported grids names as "name1, name2, name3"
  std::string print_supported_grids () const {
    const auto& grids = get_repo ();

    std::string str;
    if (grids.size()==0) {
      return str;
    }
    auto next = ++grids.begin();
    for (auto it=grids.begin(); next!=grids.end(); ++it, ++next) {
      str += it->first + ", ";
    }
    str += next->first;
    return str;
  }

  virtual void build_grid (const std::string& grid_name) = 0;
  virtual const repo_type& get_repo () const = 0;
};

inline GridsManager::grid_ptr_type
GridsManager::get_grid(const std::string& name) const
{
  scream_require_msg (supports_grid(name),
                      "Error! Grids manager '" + this->name() + "' does not provide grid '" + name + "'.\n"
                      "       Supported grids are: " + print_supported_grids()  + "\n");

  return get_repo().at(name);
}

// Forward declarations
class Comm;
class ParameterList;

// A short name for the factory for grid managers
using GridsManagerFactory 
    = util::Factory<GridsManager,
                    util::CaseInsensitiveString,
                    std::shared_ptr<GridsManager>,
                    const Comm&,const ParameterList&>;

} // namespace scream

#endif // SCREAM_GRIDS_MANAGER_HPP
