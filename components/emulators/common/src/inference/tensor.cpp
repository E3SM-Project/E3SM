/**
 * @file tensor.cpp
 * @brief Tensor and TensorMap implementation.
 */

#include "tensor.hpp"

#include "inference_error.hpp"

#include <utility>

namespace emulator {
namespace inference {

namespace {

std::int64_t product(const std::vector<std::int64_t> &dims) {
  std::int64_t n = 1;
  for (std::int64_t d : dims) {
    EMULATOR_INFER_REQUIRE(
        d >= 0, "Tensor dimensions must be non-negative, got " << d << ".");
    n *= d;
  }
  return n;
}

} // namespace

Tensor::Tensor(std::string name, std::vector<std::int64_t> dims)
    : m_name(std::move(name)), m_dims(std::move(dims)),
      m_size(product(m_dims)), m_writable(true) {
  m_storage.assign(static_cast<std::size_t>(m_size), 0.0);
  m_data = m_storage.data();
}

Tensor Tensor::make_view(std::string name, const double *data,
                         std::vector<std::int64_t> dims, bool writable) {
  Tensor t;
  t.m_name = std::move(name);
  t.m_dims = std::move(dims);
  t.m_size = product(t.m_dims);
  EMULATOR_INFER_REQUIRE(data != nullptr || t.m_size == 0,
                         "Null pointer for non-empty tensor view '" << t.m_name
                                                                    << "'.");
  t.m_data = data;
  t.m_writable = writable;
  return t;
}

Tensor Tensor::view(std::string name, double *data,
                    std::vector<std::int64_t> dims) {
  return make_view(std::move(name), data, std::move(dims), true);
}

Tensor Tensor::const_view(std::string name, const double *data,
                          std::vector<std::int64_t> dims) {
  return make_view(std::move(name), data, std::move(dims), false);
}

double *Tensor::data() {
  EMULATOR_INFER_REQUIRE(m_writable,
                         "Tensor '" << m_name << "' is a read-only view.");
  // Writable tensors are only ever built from non-const memory.
  return const_cast<double *>(m_data);
}

std::string Tensor::to_string() const {
  std::string s = m_name + "[";
  for (std::size_t i = 0; i < m_dims.size(); ++i) {
    s += (i > 0 ? "," : "") + std::to_string(m_dims[i]);
  }
  return s + "]";
}

void TensorMap::add(Tensor tensor) {
  for (const auto &t : m_tensors) {
    EMULATOR_INFER_REQUIRE(t.name() != tensor.name(),
                           "Duplicate tensor name '" << tensor.name() << "'.");
  }
  m_tensors.push_back(std::move(tensor));
}

void TensorMap::wrap(const std::string &name, double *data,
                     std::vector<std::int64_t> dims) {
  add(Tensor::view(name, data, std::move(dims)));
}

void TensorMap::wrap(const std::string &name, const double *data,
                     std::vector<std::int64_t> dims) {
  add(Tensor::const_view(name, data, std::move(dims)));
}

} // namespace inference
} // namespace emulator
