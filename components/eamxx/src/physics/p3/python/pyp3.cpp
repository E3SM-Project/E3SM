#include "physics/p3/eamxx_p3_process_interface.hpp"
#include "pygrid.hpp"
#include "pyfield.hpp"
#include "pyatmproc.hpp"

#include <pybind11/pybind11.h>

namespace scream {

namespace py = pybind11;

struct PyP3 : public PyAtmProc {
  PyP3 (const PyGrid& pygrid)
   : PyAtmProc(pygrid)
  {
    ekat::Comm comm(MPI_COMM_WORLD);
    ekat::ParameterList pl;
    pl.set("max_total_ni",10.0);

    ap = create_atmosphere_process<P3Microphysics>(comm,pl);

    // Create a grids manager on the fly
    auto gm = std::make_shared<SingleGridGM>(pygrid.grid);
    ap->set_grids(gm);
  }
};

PYBIND11_MODULE (pyp3,m) {

  m.doc() = "Python interfaces to P3Microphysics";

  py::class_<PyP3,PyAtmProc>(m,"P3")
    .def(py::init<const PyGrid&>());

}

} // namespace scream
