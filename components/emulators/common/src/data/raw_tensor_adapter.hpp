/**
 * @file raw_tensor_adapter.hpp
 * @brief Concrete TensorAdapter that hands off raw double* buffers.
 *
 * RawTensorAdapter is the always-available baseline adapter. It treats
 * the DataView buffer as-is (no data movement for to_tensor/from_tensor
 * when the DataView layout already matches), and optionally supports a
 * user-specified reshape for 4D CNN inputs or other custom shapes.
 */

#ifndef E3SM_EMULATOR_DATA_RAW_TENSOR_ADAPTER_HPP
#define E3SM_EMULATOR_DATA_RAW_TENSOR_ADAPTER_HPP

#include "tensor_adapter.hpp"

#include <cstdint>
#include <vector>

namespace emulator {

/**
 * @brief Concrete TensorAdapter using raw double* buffers.
 *
 * ## Usage
 * ```cpp
 * RawTensorAdapter adapter;
 * auto shape = adapter.tensor_shape(view);
 * const double* ptr = adapter.data_ptr(view);   // zero-copy
 *
 * // Or with a custom shape for 4D CNN:
 * adapter.set_custom_shape({1, nfields, ny, nx});
 * ```
 *
 * When a custom shape is set, tensor_shape() returns that shape and
 * data_ptr() still returns the raw pointer (reinterpretation only,
 * no data movement), provided the total element count matches.
 */
class RawTensorAdapter : public TensorAdapter {
public:
  RawTensorAdapter() = default;
  ~RawTensorAdapter() override = default;

  /**
   * @brief Set a custom tensor shape.
   *
   * The product of all dimensions must equal DataView::size() when
   * tensor_shape() / to_tensor() / data_ptr() are called.
   *
   * @param shape Desired tensor dimensions (e.g. {1, C, H, W})
   */
  void set_custom_shape(std::vector<int64_t> shape);

  /** @brief Clear any custom shape; revert to default 2D shape. */
  void clear_custom_shape();

  /** @brief Whether a custom shape is active. */
  bool has_custom_shape() const { return !m_custom_shape.empty(); }

  // ── TensorAdapter interface ─────────────────────────────────────

  /**
   * @brief Returns the tensor shape.
   *
   * Default: {dim0, dim1} from DataView::shape().
   * Custom: the shape set via set_custom_shape().
   */
  std::vector<int64_t> tensor_shape(const DataView &view) const override;

  /**
   * @brief Copy DataView data into out_buffer (a double*).
   *
   * Performs a flat memcpy of DataView::size() doubles.
   */
  void to_tensor(const DataView &view, void *out_buffer) const override;

  /**
   * @brief Copy data from in_buffer (a const double*) into the DataView.
   *
   * Performs a flat memcpy of DataView::size() doubles.
   */
  void from_tensor(const void *in_buffer, DataView &view) const override;

  /**
   * @brief Zero-copy access — always returns DataView::data().
   *
   * When a custom shape is set, validates that the total element count
   * matches. Returns nullptr on mismatch.
   */
  const double *data_ptr(const DataView &view) const override;

private:
  std::vector<int64_t> m_custom_shape;
};

} // namespace emulator

#endif // E3SM_EMULATOR_DATA_RAW_TENSOR_ADAPTER_HPP
