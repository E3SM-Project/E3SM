/**
 * @file raw_tensor_adapter.cpp
 * @brief Implementation of RawTensorAdapter.
 */

#include "raw_tensor_adapter.hpp"

#include <cstring>
#include <numeric>
#include <stdexcept>

namespace emulator {

void RawTensorAdapter::set_custom_shape(std::vector<int64_t> shape) {
  if (shape.empty()) {
    throw std::runtime_error(
        "RawTensorAdapter::set_custom_shape(): shape is empty");
  }
  m_custom_shape = std::move(shape);
}

void RawTensorAdapter::clear_custom_shape() {
  m_custom_shape.clear();
}

std::vector<int64_t>
RawTensorAdapter::tensor_shape(const DataView &view) const {
  if (!m_custom_shape.empty()) {
    // Validate total size matches
    int64_t total = std::accumulate(m_custom_shape.begin(),
                                    m_custom_shape.end(),
                                    int64_t{1},
                                    std::multiplies<int64_t>());
    if (total != static_cast<int64_t>(view.size())) {
      throw std::runtime_error(
          "RawTensorAdapter::tensor_shape(): custom shape total (" +
          std::to_string(total) + ") != DataView size (" +
          std::to_string(view.size()) + ")");
    }
    return m_custom_shape;
  }
  auto [d0, d1] = view.shape();
  return {static_cast<int64_t>(d0), static_cast<int64_t>(d1)};
}

void RawTensorAdapter::to_tensor(const DataView &view,
                                 void *out_buffer) const {
  if (!out_buffer) {
    throw std::runtime_error(
        "RawTensorAdapter::to_tensor(): null output buffer");
  }
  std::memcpy(out_buffer, view.data(), view.size() * sizeof(double));
}

void RawTensorAdapter::from_tensor(const void *in_buffer,
                                   DataView &view) const {
  if (!in_buffer) {
    throw std::runtime_error(
        "RawTensorAdapter::from_tensor(): null input buffer");
  }
  std::memcpy(view.data(), in_buffer, view.size() * sizeof(double));
}

const double *RawTensorAdapter::data_ptr(const DataView &view) const {
  if (!m_custom_shape.empty()) {
    int64_t total = std::accumulate(m_custom_shape.begin(),
                                    m_custom_shape.end(),
                                    int64_t{1},
                                    std::multiplies<int64_t>());
    if (total != static_cast<int64_t>(view.size())) {
      return nullptr; // shape mismatch → cannot do zero-copy
    }
  }
  return view.data();
}

} // namespace emulator
