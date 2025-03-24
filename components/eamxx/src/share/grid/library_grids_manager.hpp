#ifndef EAMXX_LIBRARY_GRIDS_MANAGER_HPP
#define EAMXX_LIBRARY_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

namespace scream {

// This class is meant to be used within scopes that need a grids manager
// object, but they only have pre-built grids. The user can then simply
// create a LibraryGridsManager, and add the pre-existing grids. Afterwards,
// it's business as usual with GridsManager's interfaces.

class LibraryGridsManager : public GridsManager
{
public:
  template<typename... Pointers>
  explicit LibraryGridsManager(Pointers&&... ptrs) {
    add_grids(std::forward<Pointers>(ptrs)...);
  }

  virtual ~LibraryGridsManager () = default;

  std::string name () const { return "Library grids_manager"; }

  void build_grids () override {}

  void add_grids () {}

  template<typename... Pointers>
  void add_grids (grid_ptr_type p, Pointers&&... ptrs) {
    add_grid(p);
    add_grids(std::forward<Pointers>(ptrs)...);
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
    return nullptr;
  }
};

} // namespace scream

#endif // EAMXX_LIBRARY_GRIDS_MANAGER_HPP
