/**
 * @file tensor_adapter.hpp
 * @brief Abstract interface for reshaping DataView data into tensors.
 *
 * TensorAdapter bridges DataView's flat or 2D buffer with the specific
 * tensor shapes and memory formats expected by ML/AI backends (LibTorch,
 * ONNX Runtime, raw pointers, etc.).
 *
 * The void* interface provides type erasure across backends. Each concrete
 * adapter casts internally to its native tensor type.
 */

#ifndef E3SM_EMULATOR_DATA_TENSOR_ADAPTER_HPP
#define E3SM_EMULATOR_DATA_TENSOR_ADAPTER_HPP

#include "data_view.hpp"

#include <cstdint>
#include <vector>

namespace emulator {

/**
 * @brief Abstract interface for converting DataView ↔ tensor.
 *
 * ## Design
 * - `to_tensor()` copies/reshapes DataView data into a backend tensor.
 * - `from_tensor()` copies tensor data back into the DataView.
 * - `data_ptr()` provides zero-copy access when the DataView layout
 *   already matches the desired tensor layout (returns nullptr otherwise).
 * - `tensor_shape()` reports the dimensions the adapter will produce.
 *
 * ## Backend registration
 * Concrete adapters (e.g. LibTorchAdapter, OnnxAdapter) live in separate
 * optional libraries built when the corresponding dependency is found.
 * RawTensorAdapter is always available and treats the buffer as-is.
 */
class TensorAdapter {
public:
  virtual ~TensorAdapter() = default;

  /**
   * @brief Report the tensor dimensions this adapter will produce.
   *
   * @param view The source DataView
   * @return Tensor shape as a vector of dimension sizes.
   *         e.g. {npoints, nfields} or {1, nfields, ny, nx} for 4D CNN.
   */
  virtual std::vector<int64_t>
  tensor_shape(const DataView &view) const = 0;

  /**
   * @brief Copy / reshape DataView data into a tensor buffer.
   *
   * @param view       Source DataView
   * @param out_buffer Destination — interpretation is backend-specific.
   *                   For RawTensorAdapter this is a double*.
   *                   For LibTorch this would be an at::Tensor*.
   */
  virtual void to_tensor(const DataView &view,
                         void *out_buffer) const = 0;

  /**
   * @brief Copy tensor data back into the DataView.
   *
   * @param in_buffer Source tensor — backend-specific type via void*.
   * @param view      Destination DataView (modified in-place).
   */
  virtual void from_tensor(const void *in_buffer,
                           DataView &view) const = 0;

  /**
   * @brief Zero-copy access when the DataView layout matches the tensor.
   *
   * @param view Source DataView
   * @return Pointer to the DataView's raw buffer if no reshaping is
   *         needed, or nullptr if the layout requires transformation.
   */
  virtual const double *data_ptr(const DataView &view) const = 0;
};

} // namespace emulator

#endif // E3SM_EMULATOR_DATA_TENSOR_ADAPTER_HPP
