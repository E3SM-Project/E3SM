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

    // Get array shape and strides.
    // NOTE: since the field may be padded, the strides do not necessarily
    //       match the dims. Also, the strides must be grabbed from the
    //       actual view, since the layout doesn't know them.
    pybind11::array::ShapeContainer shape (fid.get_layout().dims());
    std::vector<ssize_t> strides;

    pybind11::dtype dt;
    switch (fid.data_type()) {
      case DataType::IntType:
        dt = get_dt_and_set_strides<int>(strides);
        break;
      case DataType::FloatType:
        dt = get_dt_and_set_strides<float>(strides);
        break;
      case DataType::DoubleType:
        dt = get_dt_and_set_strides<double>(strides);
        break;
      default:
        EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
    }

    // NOTE: you MUST set the parent handle, or else you won't have view semantic
    auto data = f.get_internal_view_data_unsafe<void,Host>();
    auto this_obj = pybind11::cast(this);
    return pybind11::array(dt,shape,strides,data,pybind11::handle(this_obj));
  }

  void sync_to_host () {
    f.sync_to_host();
  }
  void sync_to_dev () {
    f.sync_to_dev();
  }
  void print() const {
    print_field_hyperslab(f);
  }
private:

  template<typename T>
  pybind11::dtype get_dt_and_set_strides (std::vector<ssize_t>& strides) const
  {
    strides.resize(f.rank());
    switch (f.rank()) {
      case 1:
      {
        auto v = f.get_view<const T*,Host>();
        strides[0] = v.stride(0)*sizeof(T);
        break;
      }
      case 2:
      {
        auto v = f.get_view<const T**,Host>();
        strides[0] = v.stride(0)*sizeof(T);
        strides[1] = v.stride(1)*sizeof(T);
        break;
      }
      case 3:
      {
        auto v = f.get_view<const T***,Host>();
        strides[0] = v.stride(0)*sizeof(T);
        strides[1] = v.stride(1)*sizeof(T);
        strides[2] = v.stride(2)*sizeof(T);
        break;
      }
      case 4:
      {
        auto v = f.get_view<const T****,Host>();
        strides[0] = v.stride(0)*sizeof(T);
        strides[1] = v.stride(1)*sizeof(T);
        strides[2] = v.stride(2)*sizeof(T);
        strides[3] = v.stride(3)*sizeof(T);
        break;
      }
      case 5:
      {
        auto v = f.get_view<const T*****,Host>();
        strides[0] = v.stride(0)*sizeof(T);
        strides[1] = v.stride(1)*sizeof(T);
        strides[2] = v.stride(2)*sizeof(T);
        strides[3] = v.stride(3)*sizeof(T);
        strides[4] = v.stride(4)*sizeof(T);
        break;
      }
      default:
        EKAT_ERROR_MSG (
            "Unsupported field rank in PyField.\n"
            " - field name: " + f.name() + "\n"
            " - field rnak: " + std::to_string(f.rank()) + "\n");
    }

    return pybind11::dtype::of<T>();
  }
};

inline void pybind_pyfield (pybind11::module& m) {
  // Field class
  pybind11::class_<PyField>(m,"Field")
    .def(pybind11::init<>())
    .def("get",&PyField::get)
    .def("sync_to_host",&PyField::sync_to_host)
    .def("sync_to_dev",&PyField::sync_to_dev)
    .def("print",&PyField::print);
}

} // namespace scream

#endif // PYFIELD_HPP
