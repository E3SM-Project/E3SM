#ifndef EAMXX_PY_FIELD_HPP
#define EAMXX_PY_FIELD_HPP

#include "share/field/field.hpp"
#include "share/core/eamxx_pysession.hpp"

#include <ekat_assert.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifdef EAMXX_HAS_TORCH
#include <torch/torch.h>
#endif

namespace scream
{

template<typename Strides, typename T,HostOrDevice HD>
void get_ptr_and_strides (Strides& strides, const void*& p, const Field& f)
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

#ifdef EAMXX_HAS_PYTORCH

torch::ScalarType get_torch_dt(const DataType dt)
{
  switch (dt) {
    case DataType::IntType:
      return torch::kInt32;
    case DataType::FloatType:
      return torch::kFloat32;
    case DataType::DoubleType:
      return torch::kFloat64;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
  }
}

template<HostOrDevice HD>
torch::Tensor create_torch_tensor (const Field& f)
{
  auto t_dev = torch::kCPU;
  int dev_idx = -1;
  if (HD==Device and ekat::OnGpu<DefaultDevice::execution_space>::value) {
#ifdef KOKKOS_ENABLE_CUDA
    t_dev = torch::kCUDA;
#elif defined(KOKKOS_ENABLE_HIP)
    t_dev = torch::kHIP;
#else
#error "Unrecognized/unsupported device for torch interoperability"
#endif
    dev_idx = DefaultDevice::execution_space().impl_idx();
  }
  torch::ScalarType dt;

  const auto& fh  = f.get_header();
  const auto& fid = fh.get_identifier();
  for (auto d : fid.get_layout().dims()) {
    shape.push_back(d);
  }

  // TODO: must check if strides are bytes or number of entries
  std::vector<int> strides
  switch (fid.data_type()) {
    case DataType::IntType:
      dt = torch::kInt32;
      get_ptr_and_strides<int,HD>(strides,data,f);
      break;
    case DataType::FloatType:
      dt = torch::kFloat32;
      get_ptr_and_strides<float,HD>(strides,data,f);
      break;
    case DataType::DoubleType:
      dt = torch::kFloat64;
      get_ptr_and_strides<double,HD>(strides,data,f);
      break;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
  }

  auto options = torch::TensorOptions()
    .dtype(get_torch_dt(f.data_type()))
    .layout(torch::kStrided)
    .device(t_dev, dev_ids);

  auto tensor = torch::from_blob(data,sizes,strides,options);
}
#endif

} // namespace scream

#endif // EAMXX_PY_FIELD_HPP
