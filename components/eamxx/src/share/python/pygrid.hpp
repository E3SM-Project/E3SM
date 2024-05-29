#ifndef PYGRID_HPP
#define PYGRID_HPP

#include "share/grid/grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "pyfield.hpp"

#include <pybind11/pybind11.h>

namespace scream {

// Small grids manager class, to hold a pre-built grid
// We will use this to build a GM on the fly from a single grid
class SingleGridGM : public GridsManager
{
public:
  SingleGridGM (const std::shared_ptr<AbstractGrid>& grid)
  {
    add_grid(grid);
  }

  std::string name () const override { return "SingleGridGM"; }

  void build_grids () override {}

protected:
  remapper_ptr_type
  do_create_remapper (const grid_ptr_type /* from_grid */,
                      const grid_ptr_type /* to_grid */) const
  {
    EKAT_ERROR_MSG ("Error! do_create_remapper not implemented for SingleGridGM.\n");
  }
};

struct PyGrid {
  std::shared_ptr<AbstractGrid> grid;

  PyGrid(const std::string& name, int ncols, int nlevs)
  {
    ekat::Comm comm(MPI_COMM_WORLD);
    grid = create_point_grid(name,ncols,nlevs,comm);
  }

  // 2d scalar, managed/unmanaged



};

} // namespace scream

#endif // PYGRID_HPP
