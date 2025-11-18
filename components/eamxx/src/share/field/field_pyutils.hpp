#ifndef EAMXX_PY_FIELD_HPP
#define EAMXX_PY_FIELD_HPP

#include "share/field/field.hpp"
#include "share/core/eamxx_pysession.hpp"

#include <ekat_assert.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <dlpack/dlpack.h>

namespace scream
{

template<typename T,HostOrDevice HD>
void get_ptr_and_strides (std::vector<ssize_t>& strides, const void*& p, const Field& f)
{
  const auto& fl  = f.get_header().get_identifier().get_layout();
  strides.resize(fl.rank());
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.get_view<const T*,HD>();
        strides[0] = v.stride_0()*sizeof(T);
        p = v.data();
        break;
      }
    case 2:
      {
        auto v = f.get_view<const T**,HD>();
        strides[0] = v.stride_0()*sizeof(T);
        strides[1] = v.stride_1()*sizeof(T);
        p = v.data();
        break;
      }
      break;
    case 3:
      {
        auto v = f.get_view<const T***,HD>();
        strides[0] = v.stride_0()*sizeof(T);
        strides[1] = v.stride_1()*sizeof(T);
        strides[2] = v.stride_2()*sizeof(T);
        p = v.data();
        break;
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank " + std::to_string(fl.rank()) + " not supported.\n");
  }
}

template<HostOrDevice HD>
pybind11::array create_py_field (const Field& f)
{
  namespace py = pybind11;

  EKAT_REQUIRE_MSG (PySession::get().is_initialized(),
      "Error! You have not initialized the Python session yet.\n");

  // Ensure numpy has been initialized
  auto numpy = PySession::get().safe_import("numpy");

  const auto& fh  = f.get_header();
  const auto& fid = fh.get_identifier();

  // Get array shape and strides (py requires size_t for shape and int64_t for strides)
  // NOTE: since the field may be padded or a subfield, so strides do not necessarily
  //       match the field dims. Also, the strides must be grabbed from the
  //       actual view, since the layout doesn't know them.
  std::vector<ssize_t> strides;
  std::vector<ssize_t>  shape;
  for (auto d : fid.get_layout().dims()) {
    shape.push_back(d);
  }

  // Also the ptr must be grabbed from the view, to make sure it points to the 1st
  // element of the field. In case f is a subfield, this is not the same as the
  // internal data pointer of the field.
  py::dtype dt;
  const void* data;
  switch (fid.data_type()) {
    case DataType::IntType:
      dt = py::dtype::of<int>();
      get_ptr_and_strides<int,HD>(strides,data,f);
      break;
    case DataType::FloatType:
      dt = py::dtype::of<float>();
      get_ptr_and_strides<float,HD>(strides,data,f);
      break;
    case DataType::DoubleType:
      dt = py::dtype::of<double>();
      get_ptr_and_strides<double,HD>(strides,data,f);
      break;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
  }

  return py::array(dt, shape, strides, data, py::none());
}

std::vector<std::int64_t> get_strides (const FieldHeader& fh)
{
  const auto& fid = fh.get_identifier();
  const auto& fl  = fid.get_layout();
  const int rank = fl.rank();

  std::vector<std::int64_t> strides(rank);
  if (rank==0) {
    return strides;
  }

  const auto& fap = fh.get_alloc_properties ();
  auto p = fh.get_parent();
  if (p) {
    strides = get_strides(*p);
    const auto& sv_info = fap.get_subview_info();
    strides.erase(strides.begin()+sv_info.dim_idx);
  } else {
    auto dims = fl.dims();
    dims.back() = fap.get_last_extent();
    strides.back() = 1;
    for (int i=rank-1; i>=1; --i) {
      strides[i-1] = strides[i]*dims[i];
    }
  }

  return strides;
}

// Get offset of field actual data from the start of the internal data pointer
std::uint64_t get_offset (const FieldHeader& fh)
{
  auto p = fh.get_parent();
  if (not p) {
    return 0;
  }

  std::uint64_t p_offset = get_offset(*p);

  const auto& fap = fh.get_alloc_properties ();
  const auto& sv_info = fap.get_subview_info();
  if (sv_info.dim_idx>0) {
    return p_offset;
  }

  return p_offset*get_strides(fh)[0]*sv_info.slice_idx;
}

// NOTE: the field is non-const, as we may need to set some extra data in the header
template<HostOrDevice HD>
DLTensor create_dl_tensor (Field& f)
{
  EKAT_REQUIRE_MSG (PySession::get().is_initialized(),
      "Error! You have not initialized the Python session yet.\n");

  DLTensor tensor;
  tensor.ndim = f.rank();

  // FieldHeader is non-const, as we may add extra data
        auto& fh  = f.get_header();
  const auto& fid = fh.get_identifier();
  const auto& fl  = fid.get_layout();

  if (not fh.has_extra_data("dlpack_shape")) {
    std::vector<std::int64_t> shape;
    for (auto d : fl.dims()) {
      shape.push_back(d);
    }
    fh.set_extra_data("dlpack_shape",shape);
  }
  if (not f.get_header().has_extra_data("dlpack_strides")) {
    std::vector<std::int64_t> strides = get_strides(fh);
    fh.set_extra_data("dlpack_strides",strides);
  }

  tensor.shape   = fh.get_extra_data<std::vector<std::int64_t>>("dlpack_shape").data();
  tensor.strides = fh.get_extra_data<std::vector<std::int64_t>>("dlpack_strides").data();
  tensor.byte_offset = get_offset(fh);

  // Also the ptr must be grabbed from the view, to make sure it points to the 1st
  // element of the field. In case f is a subfield, this is not the same as the
  // internal data pointer of the field.
  switch (fid.data_type()) {
    case DataType::IntType:
      tensor.dtype.code = kDLInt;
      tensor.dtype.bits = 8*sizeof(int);
      tensor.data = f.get_internal_view_data<int,HD>();
      break;
    case DataType::FloatType:
      tensor.dtype.code = kDLFloat;
      tensor.dtype.bits = 8*sizeof(float);
      tensor.data = f.get_internal_view_data<float,HD>();
      break;
    case DataType::DoubleType:
      tensor.dtype.code = kDLFloat;
      tensor.dtype.bits = 8*sizeof(double);
      tensor.data = f.get_internal_view_data<double,HD>();
      break;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
  }

  return tensor;
}

} // namespace scream

#endif // EAMXX_PY_FIELD_HPP
