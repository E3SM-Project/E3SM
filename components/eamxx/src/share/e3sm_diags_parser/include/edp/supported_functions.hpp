#ifndef EDP_SUPPORTED_FUNCTIONS_HPP
#define EDP_SUPPORTED_FUNCTIONS_HPP

#include <array>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
namespace edp {
struct SupportedFunction {
  std::string_view name;
  std::string_view desc;
  std::span<const std::string_view> arguments;

  std::string to_string() const {
    std::string str_{name};
    str_ += "(";
    if (!arguments.empty()) {
      for (auto val : arguments) {
        str_ += std::string(val) + ",";
      }
      str_.pop_back(); // removes trailing comma
    }
    str_ += ")";
    str_ += "\n--- " + std::string(desc);

    return str_;
  }
};
inline std::ostream& operator<<(std::ostream& os,
                                const SupportedFunction& function) {
  return os << function.to_string() << '\n';
}

inline constexpr std::array where_args{
    std::string_view{"<boolean expression>"},
};

inline constexpr std::array sum_args{
    std::string_view{"dims=[..]"},
};

inline constexpr std::array derivative_args{
    std::string_view{"dx"},
    std::string_view{"dims=[..]"},
};

inline constexpr std::array<std::string_view, 0> tend_args{};

inline constexpr std::array supported{
    SupportedFunction{
        .name = "where",
        .desc = "applies condition to operand",
        .arguments = where_args,
    },
    SupportedFunction{
        .name = "sum",
        .desc = "sums operand over designated indices (int or name)",
        .arguments = sum_args,
    },
    SupportedFunction{
        .name = "derivative",
        .desc = "takes derivative w.r.t. `dx` over designated dimension",
        .arguments = derivative_args,
    },
    SupportedFunction{
        .name = "tend",
        .desc = "calculates the tendency of a variable over time",
        .arguments = tend_args,
    },
};
} // namespace edp

#endif
