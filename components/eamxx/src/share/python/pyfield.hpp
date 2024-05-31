#ifndef PYFIELD_HPP
#define PYFIELD_HPP

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace scream {

struct PyField {
  Field f;

  // Create empty field
  PyField () = default;

  PyField(const FieldIdentifier& fid,
          const int pack_size = 1)
  {
    f = Field(fid);
    f.get_header().get_alloc_properties().request_allocation(pack_size);
    f.allocate_view();
  }

  pybind11::array get () const {

    const auto& fh  = f.get_header();
    const auto& fid = fh.get_identifier();

    // Can this actually happen? For now, no, since we only create fields from identifiers, so each PyField
    // holds separate memory. However, this may change if we allow subfields.
    EKAT_REQUIRE_MSG (f.get_header().get_parent().lock()==nullptr,
        "Error! Cannot get the array for a field that is a subfield of another. Please, get array of parent field.\n"
        "  - field name : " + fid.name() + "\n"
        "  - parent name: " + fh.get_parent().lock()->get_identifier().name() + "\n");
    pybind11::dtype dt;

    switch (fid.data_type()) {
      case DataType::IntType:     dt = pybind11::dtype::of<int>();    break;
      case DataType::FloatType:   dt = pybind11::dtype::of<float>();  break;
      case DataType::DoubleType:  dt = pybind11::dtype::of<double>(); break;
      default:                    EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
    }

    pybind11::array::ShapeContainer shape (fid.get_layout().dims());
    auto data = f.get_internal_view_data_unsafe<void,Host>();
    // auto data = reinterpret_cast<void*>(f.get_internal_view_data_unsafe<char,Host>());
    return pybind11::array(dt,shape,data);
  }

  void sync_to_host () {
    f.sync_to_host();
  }
  void print() const {
    print_field_hyperslab(f);
  }
};

inline void pybind_pyfield (pybind11::module& m) {
  // Field class
  pybind11::class_<PyField>(m,"Field")
    .def(pybind11::init<>())
    .def("get",&PyField::get)
    .def("sync_to_host",&PyField::sync_to_host)
    .def("print",&PyField::print);
}

} // namespace scream

#endif // PYFIELD_HPP
