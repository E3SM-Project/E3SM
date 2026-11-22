#ifndef E3SM_EMULATOR_FIELD_TYPES_HPP
#define E3SM_EMULATOR_FIELD_TYPES_HPP

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace e3sm::emulator {

inline constexpr std::string_view FLOAT64 = "float64";
inline constexpr std::string_view INT32 = "int32";

enum class ScalarType : int {
  Real64,
  Int32,
};

struct TypeEntry {
  std::string_view name;
  ScalarType type;
};

inline constexpr std::array dtype_table{TypeEntry{FLOAT64, ScalarType::Real64},
                                        TypeEntry{INT32, ScalarType::Int32}};

inline ScalarType dtype_string_to_enum(std::string_view dtype) {
  for (const TypeEntry& entry : dtype_table) {
    if (entry.name == dtype) {
      return entry.type;
    }
  }
  throw std::runtime_error("Unsupport type: " + std::string(dtype));
}

inline std::string_view dtype_enum_to_string(ScalarType type) {
  for (const TypeEntry& entry : dtype_table) {
    if (entry.type == type) {
      return entry.name;
    }
  }
  // Is this even possible? I suppose enums decay to ints
  throw std::runtime_error("Invalid ENUM value");
}

/**
 *Variadic template used to construct the variant from a sinlge list of types
 inorderto minimize the number of edits required to support a new type.
 */
template <typename... Ts> struct FieldTypes {

  // compare each of the supported types `Ts` with `T`
  // and combine with logical or.
  template <typename T>
  static constexpr bool contains =
      (std::same_as<std::remove_cv_t<T>, Ts> || ...);

  // The Variants: span applies const is IsConst=true.
  template <bool IsConst>
  using Views =
      std::variant<std::span<std::conditional_t<IsConst, const Ts, Ts>>...>;

  using Storage = std::variant<std::vector<Ts>...>;
};

// How to use SupportedTypes:
// - SupportedTypes::contains<double> checks if `double` is Supported.
// - SupportedTypes::Views<false> -> variant<span<double>,span<int32_t>>
//
using SupportedTypes = FieldTypes<double, std::int32_t>;

template <typename T> struct ScalarTypeOf;

template <> struct ScalarTypeOf<double> {
  static constexpr auto value = ScalarType::Real64;
};

template <> struct ScalarTypeOf<std::int32_t> {
  static constexpr auto value = ScalarType::Int32;
};
// This is a helper to get ScalarType from T
template <typename T>
inline constexpr ScalarType scalar_type_of =
    ScalarTypeOf<std::remove_cv_t<T>>::value;

template <typename T>
concept FieldScalar = SupportedTypes::contains<T>;

template <typename T, bool IsConst>
concept CompatibleElement = FieldScalar<T> && (IsConst || !std::is_const_v<T>);

inline SupportedTypes::Storage make_storage_vector(ScalarType type,
                                                   std::size_t size) {
  switch (type) {
  case ScalarType::Real64:
    return std::vector<double>(size);
  case ScalarType::Int32:
    return std::vector<std::int32_t>(size);
  }
  throw std::invalid_argument("Unsupported field scalar type");
}

inline SupportedTypes::Views<false>
make_field_view(ScalarType type, std::size_t size, void* ptr) {
  if (ptr == nullptr || size == 0) {
    throw std::runtime_error("Can't make span with nullptr or 0 extent");
  }
  switch (type) {
  case ScalarType::Real64:
    return std::span{static_cast<double*>(ptr), size};
  case ScalarType::Int32:
    return std::span{static_cast<std::int32_t*>(ptr), size};
  }
  throw std::invalid_argument("Unsupported field scalar type");
}
/**
 * BasicFieldView
 *  This class acts as a wrapper around a std::variant that contains
 *  all the supported types to be used for as a field for an emulator, coupling,
 * etc...
 *
 *  `using FieldView = BasicFieldView<false>` for mutable span<T>
 *  `using ConstFieldView = BasicFieldView<true>` for read-only span<const T>
 *
 *
 */
template <bool IsConst> class BasicFieldView {

  // Element = const T if IsConst is true : T otherwise
  template <typename T> using Element = std::conditional_t<IsConst, const T, T>;
  using Storage = SupportedTypes::Views<IsConst>;

public:
  // Constructor from span
  template <typename T, std::size_t Extent>
    requires CompatibleElement<T, IsConst>
  BasicFieldView(std::span<T, Extent> data)
      : data_(std::span<Element<std::remove_const_t<T>>>{data}) {}

  explicit BasicFieldView(Storage data) : data_(std::move(data)) {}

  BasicFieldView() = default;

  // Constructor from vector:
  template <typename T, typename Allocator>
    requires CompatibleElement<T, IsConst>
  BasicFieldView(std::vector<T, Allocator>& values)
      : BasicFieldView(std::span<T>{values}) {}

  // Constructor from const vector:
  template <typename T, typename Allocator>
    requires(IsConst && FieldScalar<T>)
  BasicFieldView(const std::vector<T, Allocator>& values)
      : BasicFieldView(std::span<const T>{values}) {}

  // Delete these to force named vector construction **only**
  template <typename T, typename Allocator>
  BasicFieldView(std::vector<T, Allocator>&&) = delete;
  template <typename T, typename Allocator>
  BasicFieldView(const std::vector<T, Allocator>&&) = delete;

  // Generic wrapper around std::visit
  template <typename Function> decltype(auto) visit(Function&& function) const {
    return std::visit(std::forward<Function>(function), data_);
  }

  // Helper methods that wrap visit to allow `.size` and `.empty` calls like a
  // normal span
  std::size_t size() const {
    return visit([](auto values) { return values.size(); });
  }
  bool empty() const {
    return visit([](auto values) { return values.empty(); });
  }

  // Get the pointer address of the underlying span
  const void* addr() const {
    return visit(
        [](auto values) { return static_cast<const void*>(values.data()); });
  }

  BasicFieldView subspan(std::size_t offset, std::size_t count) const {
    if (offset > size() || count > size() - offset) {
      throw std::out_of_range("FieldView::subspan: invalid range");
    }

    return visit([=](auto values) {
      return BasicFieldView{values.subspan(offset, count)};
    });
  }

  std::ostringstream print(std::size_t count) const {

    std::ostringstream os;

    if (count > 0) {

      os << "       [ ";
      visit([&](auto values) {
        const auto n = std::min(count, values.size());

        for (std::size_t i = 0; i < n; ++i) {
          if (i > 0) {
            os << ", ";
          }
          os << values[i];
        }
      });

      os << " ]\n";
      return os;
    }
  }

  // Helper to make a non-const view read-only and it only
  // exists for ConstFieldView.  For that it's signature becomes:
  // ConstFieldView(const FieldView& other)
  template <bool OtherConst>
    requires(IsConst && !OtherConst)
  BasicFieldView(const BasicFieldView<OtherConst>& other)
      : data_(other.visit([](auto values) -> Storage {
          using T = typename decltype(values)::value_type;
          return std::span<const T>{values};
        })) {}

  // Get's the corresponding ScalarType enum
  ScalarType scalar_type() const {
    return visit([](auto values) {
      using T = typename decltype(values)::value_type;
      return scalar_type_of<T>;
    });
  }

  std::size_t size_bytes() const {
    return std::visit([](auto values) { return values.size_bytes(); }, data_);
  }

  // IF YOU KNOW THE TYPE, returns the underlying span
  // We know the type in the variant during a visit call, so this can be used
  // below in the copy_from method.
  template <typename T> std::span<Element<T>> as() const {
    return std::get<std::span<Element<T>>>(data_);
  }

  bool is_same_view_as(const BasicFieldView<true>& other) const {
    return scalar_type() == other.scalar_type() && size() == other.size() &&
           addr() == other.addr();
  }

  // usage: dest.copy_from(source)
  bool copy_from(const BasicFieldView<true>& source) const
    requires(!IsConst)
  {
    if (scalar_type() != source.scalar_type()) {
      throw std::invalid_argument("Field copy: scalar type mismatch");
    }
    if (size() != source.size()) {
      throw std::invalid_argument("Field copy: size mismatch");
    }
    if (is_same_view_as(source)) {
      return false;
    }

    visit([&](auto destination) {
      using T = typename decltype(destination)::value_type;
      std::ranges::copy(source.template as<T>(), destination.begin());
    });
    return true;
  }

private:
  Storage data_;
};

using FieldView = BasicFieldView<false>;
using ConstFieldView = BasicFieldView<true>;

/**
 * FieldStorage --
 * the owning version of BasicFieldView
 */

class FieldStorage {
  using Storage = SupportedTypes::Storage;

public:
  FieldStorage() = default;
  explicit FieldStorage(ScalarType type, std::size_t size = 0)
      : data_(make_storage_vector(type, size)) {}

  // Create views from storage:
  FieldView view() {
    return std::visit([](auto& values) { return FieldView{std::span{values}}; },
                      data_);
  }
  ConstFieldView view() const {
    return std::visit(
        [](const auto& values) { return ConstFieldView{std::span{values}}; },
        data_);
  }

  std::size_t size() const {
    return std::visit([](auto values) { return values.size(); }, data_);
  }

  ScalarType scalar_type() const {
    return std::visit(
        [](const auto& values) {
          using T = typename std::remove_cvref_t<decltype(values)>::value_type;
          return scalar_type_of<T>;
        },
        data_);
  }

  void resize(std::size_t size) {
    std::visit([size](auto& values) { values.resize(size); }, data_);
  }

  void zero() {
    view().visit([](auto values) {
      using T = typename decltype(values)::value_type;
      std::ranges::fill(values, T{0});
    });
  }

  static FieldStorage like(ConstFieldView source) {
    return FieldStorage{source.scalar_type(), source.size()};
  }

private:
  Storage data_;
};

// Helper to be used when iterating over two views of the same type
template <typename Function>
void visit_matching(FieldView destination, ConstFieldView source,
                    Function&& function) {
  if (destination.scalar_type() != source.scalar_type()) {
    throw std::invalid_argument("Field scalar type mismatch");
  }
  if (destination.size() != source.size()) {
    throw std::invalid_argument("Field size mismatch");
  }

  destination.visit([&](auto dest) {
    using T = typename decltype(dest)::value_type;
    std::invoke(std::forward<Function>(function), dest, source.as<T>());
  });
}

} // namespace e3sm::emulator

#endif
