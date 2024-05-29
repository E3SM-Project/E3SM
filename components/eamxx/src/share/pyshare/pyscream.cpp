#include "scream_session.hpp"
#include "pyfield.hpp"
#include "pygrid.hpp"
#include "pyunits.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace scream {

void initialize () {
  initialize_scream_session(true);
}
void finalize () {
  finalize_scream_session();
}

PYBIND11_MODULE (pyscream,m) {

  m.doc() = "Python interfaces to certain EAMxx infrastructure code";

  // Scream Session
  m.def("init",&initialize);
  m.def("finalize",&finalize);

  // Units
  py::module units = m.def_submodule("units");
  py::class_<PyUnits>(units,"Units")
    .def(py::init<>())
    .def(py::init<const ekat::units::Units&>())
    .def("str",&PyUnits::str)
    .def("__mul__",[](const PyUnits& x, const PyUnits& y) { return PyUnits(x.units*y.units); })
    .def("__truediv__",[](const PyUnits& x, const PyUnits& y) { return PyUnits(x.units/y.units); });

  units.def("pow",[](const PyUnits& x, const int n) { return PyUnits(pow(x.units,ekat::RationalConstant(n))); });
  units.def("root",[](const PyUnits& x, const int n) { return PyUnits(pow(x.units,ekat::RationalConstant(1,n))); });
  units.attr("one") = &short_units::one;
  units.attr("m") = &short_units::m;
  units.attr("s") = &short_units::s;
  units.attr("kg") = &short_units::kg;
  units.attr("K") = &short_units::K;
  units.attr("N") = &short_units::N;
  units.attr("J") = &short_units::J;
  units.attr("W") = &short_units::W;
  units.attr("Pa") = &short_units::Pa;

  // Field class
  py::class_<PyField>(m,"Field")
    .def(py::init<>())
    .def(py::init<const std::string&,
                  const std::vector<int>&>())
    .def(py::init<const std::string&,
                  py::array_t<double>
                 >())
    .def("print",&PyField::print)
    .def("cleanup",&PyField::cleanup);

  // Grid
  py::class_<PyGrid>(m,"Grid")
    .def(py::init<>())
    .def("scalar2d",py::overload_cast<const std::string&,const PyUnits&>(&PyGrid::scalar2d))
    .def("scalar2d",py::overload_cast<const std::string&,const PyUnits&,py::array_t<double>>(&PyGrid::scalar2d))
    .def("vector2d",py::overload_cast<const std::string&,const PyUnits&,int>(&PyGrid::vector2d))
    .def("vector2d",py::overload_cast<const std::string&,const PyUnits&,int,py::array_t<double>>(&PyGrid::vector2d))
    .def("scalar3d_mid",py::overload_cast<const std::string&,const PyUnits&>(&PyGrid::scalar3d_mid))
    .def("scalar3d_mid",py::overload_cast<const std::string&,const PyUnits&,py::array_t<double>>(&PyGrid::scalar3d_mid))
    .def("scalar3d_int",py::overload_cast<const std::string&,const PyUnits&>(&PyGrid::scalar3d_int))
    .def("scalar3d_int",py::overload_cast<const std::string&,const PyUnits&,py::array_t<double>>(&PyGrid::scalar3d_int));

  // GridsManager
  py::class_<PyGridsManager>(m,"GridsManager")
    .def(py::init<const std::string&,int,int>())
    .def("get_grid",&PyGridsManager::get_grid)
    .def("cleanup",&PyGridsManager::cleanup);
}

} // namespace scream
