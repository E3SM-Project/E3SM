#ifndef EAMXX_LIBRARY_GRIDS_MANAGER_HPP
#define EAMXX_LIBRARY_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

#include "ekat/mpi/ekat_comm.hpp"

namespace scream {

// This class is meant to be used within scopes that need a grids manager
// object, but they only have 1+ pre-built grids. The user can then simply
// create a LibraryGridsManager, and add the pre-existing grids. Afterwards,
// it can rely on GridsManager's interfaces.

class LibraryGridsManager : public GridsManager
{
public:
  template<typename... GridsPointers>
  explicit LibraryGridsManager(GridPtrs&&... ptrs) {
    (add_grid(std::forward<GridPtrs>(ptrs)), ...);
  }

  virtual ~LibraryGridsManager () = default;

  std::string name () const { return "Library grids_manager"; }

  void build_grids () override {}

  // Overload, not hide
  using GridsManager::add_grid;
  void add_grid (const grid_ptr_type& grid) {
    GridsManager::add_grid (grid);
  }

protected:

  remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const
  {
    EKAT_ERROR_MSG (
        "Error! LibraryGridsManager is not capable of creating remappers.\n"
        " - from_grid: " + from_grid->name() + "\n"
        " - to_grid:   " + to_grid->name() + "\n");
        
  }
};

} // namespace scream

#endif // EAMXX_LIBRARY_GRIDS_MANAGER_HPP
