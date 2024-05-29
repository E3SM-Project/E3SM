#include "physics/p3/eamxx_p3_process_interface.hpp"
#include "pygrid.hpp"
#include "pyfield.hpp"

#include <pybind11/pybind11.h>

namespace scream {

namespace py = pybind11;

struct PyP3 {
  std::shared_ptr<AtmosphereProcess> p3;

  PyP3 (const PyGridsManager& pygm) {
    ekat::Comm comm(MPI_COMM_WORLD);
    ekat::ParameterList pl;
    pl.set("max_total_ni",10.0);

    p3 = create_atmosphere_process<P3Microphysics>(comm,pl);
    p3->set_grids(pygm.gm);
  }

  void set_fields(const std::vector<PyField>& pyfields) {
    for (const auto& pyf : pyfields) {
      const auto& fid = pyf.f.get_header().get_identifier();
      if (p3->has_required_field(fid)) {
        p3->set_required_field(pyf.f.get_const());
      }
      if (p3->has_computed_field(fid)) {
        p3->set_computed_field(pyf.f);
      }
    }
  }

  void initialize () {
    util::TimeStamp ts;
    p3->initialize(ts,RunType::Initial);
  }

  void cleanup () {
    p3 = nullptr;
  }
};

PYBIND11_MODULE (pyp3,m) {

  m.doc() = "Python interfaces to P3Microphysics";

  py::class_<PyP3>(m,"P3")
    .def(py::init<const PyGridsManager&>())
    .def("set_fields",&PyP3::set_fields)
    .def("initialize",&PyP3::initialize)
    .def("cleanup",&PyP3::cleanup);

}

} // namespace scream

