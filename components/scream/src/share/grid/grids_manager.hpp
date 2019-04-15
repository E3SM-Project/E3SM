#ifndef SCREAM_GRIDS_MANAGER_HPP
#define SCREAM_GRIDS_MANAGER_HPP

#include "share/grid/abstract_grid.hpp"
#include <set>

namespace scream
{

class GridsManager
{
public:
  using grid_type     = AbstractGrid;
  using grid_ptr_type = std::shared_ptr<grid_type>;
  using repo_type     = std::map<GridType, grid_ptr_type>;

  GridsManager () = default;
  virtual ~GridsManager () = default;

  std::shared_ptr<grid_type> get_grid (const GridType grid_type) const;

  virtual void build_grids (const std::set<GridType>& grid_types) = 0;
protected:

  virtual const repo_type& get_repo () const = 0;
};

inline GridsManager::grid_ptr_type
GridsManager::get_grid(const GridType grid_type) const
{
  auto it = get_repo().find(grid_type);
  scream_require_msg (it!=get_repo().end(), "Error! Grid type '" + e2str(grid_type) + "' not registered.\n");

  return it->second;
}

// A short name for the factory for grid managers
using GridsManagerFactory = util::Factory<GridsManager,util::CaseInsensitiveString,const ParameterList&>;

} // namespace scream

#endif // SCREAM_GRIDS_MANAGER_HPP
