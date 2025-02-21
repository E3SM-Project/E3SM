#ifndef PYFIELD_HPP
#define PYFIELD_HPP

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/list.h>

namespace nb = nanobind;

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

  template <typename FRAMEWORK>
  nb::ndarray<FRAMEWORK> get () const {
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
    int shape_t = f.rank();
    size_t shape[shape_t] = {0};
    for (int i=0; i<shape_t; ++i) {
      shape[i] = fid.get_layout().dims()[i];
    }
    std::vector<ssize_t> strides;

    nb::dlpack::dtype dt;
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
    auto this_obj = nb::cast(this);
    return nb::ndarray<FRAMEWORK>(data, shape_t, shape, nb::handle(this_obj), strides.data(), dt);
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
  nb::dlpack::dtype get_dt_and_set_strides (std::vector<ssize_t>& strides) const
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

    return nb::dtype<T>();
  }
};

inline void nb_pyfield (nb::module_& m) {
  // Field class
  nb::class_<PyField>(m,"Field")
    .def(nb::init<>())
    .def("get",&PyField::get<nb::numpy>)
    .def("sync_to_host",&PyField::sync_to_host)
    .def("sync_to_dev",&PyField::sync_to_dev)
    .def("print",&PyField::print);
}

} // namespace scream

#endif // PYFIELD_HPP
