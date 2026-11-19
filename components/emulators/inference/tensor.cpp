/**
 * @file tensor.cpp
 * @brief Tensor and TensorMap implementation.
 */

#include "tensor.hpp"
#include "inference_error.hpp"
#include <field_types.hpp>

#include <utility>

namespace e3sm::emulator::inference {

namespace {

std::int64_t product(const std::vector<std::int64_t>& dims) {
  std::int64_t n = 1;
  for (std::int64_t d : dims) {
    EMULATOR_INFER_REQUIRE(
        d >= 0, "Tensor dimensions must be non-negative, got " << d << ".");
    n *= d;
  }
  return n;
}

} // namespace

Tensor::Tensor(std::string name, std::vector<std::int64_t> dims, Buffer buffer)
    : name_(std::move(name)), dims_(std::move(dims)), size_(product(dims_)),
      data_(std::move(buffer)) {

  const auto values = std::as_const(*this).view();

  EMULATOR_INFER_REQUIRE(values.size() == static_cast<std::size_t>(size_),
                         "Tensor '" << name_ << "': shape requires " << size_
                                    << " elements, but buffer has "
                                    << values.size() << ".");
}
Tensor::Tensor(std::string name, std::vector<std::int64_t> dims,
               ScalarType type)
    : Tensor(
          std::move(name), dims,
          Buffer{FieldStorage{type, static_cast<std::size_t>(product(dims))}}) {
}

// Return the underlying View in the Buffer variant
FieldView Tensor::view() {
  return std::visit(
      [](auto& buffer) -> FieldView {
        using B = std::remove_cvref_t<decltype(buffer)>;

        if constexpr (std::same_as<B, FieldStorage>) {
          return buffer.view();
        } else if constexpr (std::same_as<B, FieldView>) {
          return buffer;
        } else {
          throw InferenceError("Tensor is a read-only view");
        }
      },
      data_);
}

ConstFieldView Tensor::view() const {
  return std::visit(
      [](const auto& buffer) -> ConstFieldView {
        using B = std::remove_cvref_t<decltype(buffer)>;

        if constexpr (std::same_as<B, FieldStorage>) {
          return buffer.view();
        } else {
          return buffer;
        }
      },
      data_);
}

std::string Tensor::to_string() const {
  std::string s = name_ + "[";
  for (std::size_t i = 0; i < dims_.size(); ++i) {
    s += (i > 0 ? "," : "") + std::to_string(dims_[i]);
  }
  return s + "]";
}

void TensorMap::add(Tensor tensor) {
  for (const auto& t : tensors_) {
    EMULATOR_INFER_REQUIRE(t.name() != tensor.name(),
                           "Duplicate tensor name '" << tensor.name() << "'.");
  }
  tensors_.push_back(std::move(tensor));
}

Tensor Tensor::from_view(std::string name, FieldView data,
                         std::vector<std::int64_t> dims) {
  return Tensor(std::move(name), dims, Buffer{data});
}

Tensor Tensor::from_const_view(std::string name, ConstFieldView data,
                               std::vector<std::int64_t> dims) {
  return Tensor(std::move(name), dims, Buffer{data});
}

void TensorMap::wrap(const std::string& name, FieldView data,
                     std::vector<std::int64_t> dims) {
  add(Tensor::from_view(name, data, std::move(dims)));
}

void TensorMap::wrap(const std::string& name, ConstFieldView data,
                     std::vector<std::int64_t> dims) {
  add(Tensor::from_const_view(name, data, std::move(dims)));
}

} // namespace e3sm::emulator::inference
