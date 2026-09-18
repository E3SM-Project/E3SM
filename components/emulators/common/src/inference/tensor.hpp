/**
 * @file tensor.hpp
 * @brief Named, shaped buffers exchanged with an inference backend.
 */

#ifndef E3SM_EMULATOR_INFERENCE_TENSOR_HPP
#define E3SM_EMULATOR_INFERENCE_TENSOR_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace emulator {
namespace inference {

/**
 * @brief A named, shaped, row-major buffer of doubles.
 *
 * A Tensor either owns its memory or views memory owned by the caller, so
 * model fields can be handed to a backend without a copy. Elements are
 * always double (E3SM's real(r8)); a backend converts to the model's
 * precision if it needs to.
 *
 * Row-major: a Fortran array a(nlev,ncol) corresponds to dims {ncol, nlev}.
 *
 * Move-only, so a large field is never copied by accident.
 */
class Tensor {
public:
  /// @brief Allocate a zero-filled tensor.
  Tensor(std::string name, std::vector<std::int64_t> dims);

  /// @brief View writable memory owned by the caller.
  static Tensor view(std::string name, double *data,
                     std::vector<std::int64_t> dims);

  /// @brief View read-only memory owned by the caller.
  static Tensor const_view(std::string name, const double *data,
                           std::vector<std::int64_t> dims);

  Tensor(const Tensor &) = delete;
  Tensor &operator=(const Tensor &) = delete;
  Tensor(Tensor &&) = default;
  Tensor &operator=(Tensor &&) = default;

  const std::string &name() const { return m_name; }
  const std::vector<std::int64_t> &dims() const { return m_dims; }

  /// @brief Total element count (product of dims).
  std::int64_t size() const { return m_size; }

  bool writable() const { return m_writable; }

  /**
   * @brief Writable pointer to the data (null if the tensor is empty).
   * @throws InferenceError if this is a read-only view
   */
  double *data();

  /// @brief Read-only pointer to the data (null if the tensor is empty).
  const double *cdata() const { return m_data; }

  /// @brief "name[d0,d1,...]", for messages.
  std::string to_string() const;

private:
  Tensor() = default;

  static Tensor make_view(std::string name, const double *data,
                          std::vector<std::int64_t> dims, bool writable);

  std::string m_name;
  std::vector<std::int64_t> m_dims;
  std::int64_t m_size = 0;
  std::vector<double> m_storage; ///< Used only when the tensor owns its data
  const double *m_data = nullptr;
  bool m_writable = false;
};

/**
 * @brief An ordered set of tensors with unique names.
 *
 * Order is preserved because backends may pass tensors to a model
 * positionally.
 */
class TensorMap {
public:
  /**
   * @brief Append a tensor.
   * @throws InferenceError if the name is already present
   */
  void add(Tensor tensor);

  /// @brief Append a view of writable caller memory.
  void wrap(const std::string &name, double *data,
            std::vector<std::int64_t> dims);

  /// @brief Append a view of read-only caller memory.
  void wrap(const std::string &name, const double *data,
            std::vector<std::int64_t> dims);

  std::size_t size() const { return m_tensors.size(); }

  std::vector<Tensor>::iterator begin() { return m_tensors.begin(); }
  std::vector<Tensor>::iterator end() { return m_tensors.end(); }
  std::vector<Tensor>::const_iterator begin() const {
    return m_tensors.begin();
  }
  std::vector<Tensor>::const_iterator end() const { return m_tensors.end(); }

private:
  std::vector<Tensor> m_tensors;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_INFERENCE_TENSOR_HPP
