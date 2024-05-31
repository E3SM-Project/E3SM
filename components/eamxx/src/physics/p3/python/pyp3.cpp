#include "physics/p3/eamxx_p3_process_interface.hpp"
#include "pygrid.hpp"
#include "pyfield.hpp"
#include "pyatmproc.hpp"
#include "pyparamlist.hpp"

#include <pybind11/pybind11.h>

namespace scream {

namespace py = pybind11;

struct PyP3 : public PyAtmProc {
  PyP3 (const PyGrid& pygrid, const PyParamList& params)
   : PyAtmProc(pygrid)
  {
    init (params.pl);
  }
  PyP3 (const PyGrid& pygrid)
   : PyAtmProc(pygrid)
  {
    // P3 does not have a default (in the C++ code) for this parameter,
    // so an entry MUST be set
    ekat::ParameterList pl;
    pl.set("max_total_ni",740.0e3);
    init(pl);
  }

  void init(const ekat::ParameterList& params) {
    auto comm = pygrid.grid->get_comm();
    ap = create_atmosphere_process<P3Microphysics>(comm,params);

    // Create a grids manager on the fly
    auto gm = std::make_shared<SingleGridGM>(pygrid.grid);
    ap->set_grids(gm);

    set_fields(create_fields());
  }
};

PYBIND11_MODULE (pyp3,m) {

  m.doc() = "Python interfaces to P3Microphysics";

  py::class_<PyP3,PyAtmProc>(m,"P3")
    .def(py::init<const PyGrid&>())
    .def(py::init<const PyGrid&,const PyParamList&>());

}

} // namespace scream
