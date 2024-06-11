#ifndef PYGRID_HPP
#define PYGRID_HPP

#include "share/grid/grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "pyfield.hpp"
#include "pyutils.hpp"

#include <pybind11/pybind11.h>

#include <mpi.h>

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
                      const grid_ptr_type /* to_grid */) const override
  {
    EKAT_ERROR_MSG ("Error! do_create_remapper not implemented for SingleGridGM.\n");
  }
};

struct PyGrid {
  std::shared_ptr<AbstractGrid> grid;

  PyGrid () = default;

  PyGrid(const std::string& name, int ncols, int nlevs)
    : PyGrid (name,ncols,nlevs,MPI_COMM_WORLD)
  {}

  PyGrid(const std::string& name, int ncols, int nlevs, pybind11::object py_comm)
    : PyGrid (name,ncols,nlevs,get_c_comm(py_comm))
  {}

  PyGrid(const std::string& name, int ncols, int nlevs, MPI_Comm mpi_comm)
  {
    ekat::Comm comm(mpi_comm);
    grid = create_point_grid(name,ncols,nlevs,comm);
  }
};

inline void pybind_pygrid (pybind11::module& m) {
  pybind11::class_<PyGrid>(m,"Grid")
    .def(pybind11::init<const std::string&,int,int>())
    .def(pybind11::init<const std::string&,int,int,pybind11::object>());
}

} // namespace scream

#endif // PYGRID_HPP
