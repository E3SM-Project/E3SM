#include "scream_session.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

using namespace scream;
using namespace ShortFieldTagsNames;
namespace py = pybind11;

void initialize () {
  initialize_scream_session(false);
}
void finalize () {
  finalize_scream_session();
}

struct PyField {
  Field f;

  // Create field and allocate memory
  PyField(const std::string& n,
          const std::vector<int>& d)
  {
    std::vector<FieldTag> t (d.size(),CMP);
    FieldIdentifier fid(n,FieldLayout(t,d),ekat::units::Units::invalid(),"");
    f = Field(fid);
    f.allocate_view();
  }

  // Create field from pre-existing python nd array
  PyField(const std::string& n,
          py::array_t<double> arr)
  {
    int rank = arr.ndim();
    std::vector<int> d(rank,-1);
    for (int n=0; n<rank; ++n) {
      d[n] = arr.shape(n);
    }
    std::vector<FieldTag> t (rank,CMP);
    FieldIdentifier fid(n,FieldLayout(t,d),ekat::units::Units::invalid(),"");
    switch (rank) {
      case 1:
      {
        Field::view_dev_t<double*> v(arr.mutable_data(0),d[0]);
        f = Field(fid,v);
        break;
      }
      case 2:
      {
        Field::view_dev_t<double**> v(arr.mutable_data(0,0),d[0],d[1]);
        f = Field(fid,v);
        break;
      }
      default:
        EKAT_ERROR_MSG ("AAARGH!\n");
    }
  }

  void print() const {
    print_field_hyperslab(f);
  }
};

PYBIND11_MODULE (pyscream,m) {

  m.doc() = "Basic interface to scream session initialization";

  m.def("init",&initialize);
  m.def("finalize",&finalize);

  py::enum_<scream::FieldTag>(m,"Tag")
    .value("INV" ,INV)
    .value("EL"  ,EL)
    .value("COL" ,COL)
    .value("GP"  ,GP)
    .value("TL"  ,TL)
    .value("LEV" ,LEV)
    .value("ILEV",ILEV)
    .value("CMP" ,CMP);

  py::class_<PyField>(m,"Field")
    .def(py::init<const std::string&,
                  const std::vector<int>&>())
    .def(py::init<const std::string&,
                  py::array_t<double>
                 >())
    .def("print",&PyField::print);
}
