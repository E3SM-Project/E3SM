#include "share/grid/point_grid.hpp"

#include <pybind11/pybind11.h>

namespace scream {

struct PyPointGrid {
  std::shared_ptr<PointGrid> g;

  // Create field and allocate memory
  PyPointGrid(const std::string& name, int ncols, int nlevs)
  {
    ekat::Comm comm(MPI_COMM_WORLD);
    g = create_point_grid (name,ncols,nlevs,comm);
  }

  void cleanup () {
    g = nullptr;
  }
};

} // namespace scream
