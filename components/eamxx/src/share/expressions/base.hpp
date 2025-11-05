#ifndef EAMXX_EXPRESSION_HPP
#define EAMXX_EXPRESSION_HPP

#include <Kokkos_Core.hpp>

#include <share/field/field.hpp>
#include <share/field/field_layout.hpp>
#include <share/core/eamxx_types.hpp>

namespace scream {

struct EvalData {
  int i = -1;
  int j = -1;
  int k = -1;
};

class ExpressionBase {
public:
  ~ExpressionBase () = default;

  virtual void evaluate (Field& result) = 0;
};

template<typename Derived>
class Expression : public ExpressionBase {
public:
  static constexpr bool is_leaf = false;

  int num_indices () const { return cast().num_indices(); }

  KOKKOS_INLINE_FUNCTION
  Real eval() const {
    return cast().eval();
  }

  void evaluate (Field& result) override;

  KOKKOS_INLINE_FUNCTION
  const Derived& cast () const { return static_cast<const Derived&>(*this); }
        Derived& cast ()       { return static_cast<      Derived&>(*this); }

  KOKKOS_INLINE_FUNCTION
  void set_eval_data (const EvalData& data) const { cast().set_eval_data(data); }
  void set_eval_layout (const FieldLayout& fl) { cast().set_eval_layout(fl); }
};

#include "evaluate.hpp"

} // namespace scream

#endif // EAMXX_EXPRESSION_HPP
