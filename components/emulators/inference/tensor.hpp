/**
 * @file tensor.hpp
 * @brief Named, shaped buffers exchanged with an inference backend.
 */

#ifndef E3SM_EMULATOR_INFERENCE_TENSOR_HPP
#define E3SM_EMULATOR_INFERENCE_TENSOR_HPP

#include <field_types.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace e3sm::emulator::inference {

using SType = e3sm::emulator::ScalarType;
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

using Buffer = std::variant<FieldStorage,  // Owned, writable storage
                            FieldView,     // Borrowed, writable storage
                            ConstFieldView // Borrowed, read-only storage
                            >;
class Tensor {

public:
  /// @brief Allocate a zero-filled owned tensor.
  Tensor(std::string name, std::vector<std::int64_t> dims,
         SType type = SType::Real64);

  /// @brief View writable memory owned by the caller.
  static Tensor from_view(std::string name, FieldView data,
                          std::vector<std::int64_t> dims);

  /// @brief View read-only memory owned by the caller.
  static Tensor from_const_view(std::string name, ConstFieldView data,
                                std::vector<std::int64_t> dims);

  Tensor(const Tensor&) = delete;
  Tensor& operator=(const Tensor&) = delete;
  Tensor(Tensor&&) = default;
  Tensor& operator=(Tensor&&) = default;

  FieldView view();
  ConstFieldView view() const;

  bool writable() const {
    return !std::holds_alternative<ConstFieldView>(data_);
  }

  ScalarType scalar_type() const { return view().scalar_type(); }

  const std::string& name() const { return name_; }
  std::span<const std::int64_t> dims() const { return dims_; }

  /// @brief Total element count (product of dims).
  std::int64_t size() const { return size_; }

  /// @brief "name[d0,d1,...]", for messages.
  std::string to_string() const;

  const void* cdata() const { return view().addr(); }
  void* data() {
    return view().visit(
        [](auto values) { return static_cast<void*>(values.data()); });
  }

  // More ergonomic access of elements for read/write:
  //  Example:  Tensor x; x.as<double>()[1] = 2.0
  template <FieldScalar T> std::span<T> as() { return view().template as<T>(); }

  template <FieldScalar T> std::span<const T> as() const {
    return view().template as<T>();
  }

private:
  // Constructor of view from external buffer
  Tensor(std::string name, std::vector<std::int64_t> dims, Buffer buffer);

  std::string name_;
  std::vector<std::int64_t> dims_;
  std::int64_t size_ = 0;
  Buffer data_;
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
  void wrap(const std::string& name, FieldView, std::vector<std::int64_t> dims);

  /// @brief Append a view of read-only caller memory.
  void wrap(const std::string& name, ConstFieldView,
            std::vector<std::int64_t> dims);

  std::size_t size() const { return tensors_.size(); }

  std::vector<Tensor>::iterator begin() { return tensors_.begin(); }
  std::vector<Tensor>::iterator end() { return tensors_.end(); }
  std::vector<Tensor>::const_iterator begin() const { return tensors_.begin(); }
  std::vector<Tensor>::const_iterator end() const { return tensors_.end(); }

private:
  std::vector<Tensor> tensors_;
};

} // namespace e3sm::emulator::inference

#endif // E3SM_EMULATOR_INFERENCE_TENSOR_HPP
