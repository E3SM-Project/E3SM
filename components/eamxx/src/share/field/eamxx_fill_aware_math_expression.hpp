#ifndef EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP
#define EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_expression_binary_op.hpp>
#include <ekat_expression_meta.hpp>

namespace ekat {

#define FILL_AWARE_BINARY_MATH_EXPR(impl,name) \
  template<typename EArg1, typename EArg2>                                                                \
  class name##FAExpression : public name##Expression<EArg1,EArg2> {                                       \
  public:                                                                                                 \
    using base_t = name##Expression<EArg1,EArg2>;                                                         \
                                                                                                          \
    static constexpr bool expr_l = is_expr_v<EArg1>;                                                      \
    static constexpr bool expr_r = is_expr_v<EArg2>;                                                      \
                                                                                                          \
    using eval_arg1_t = typename base_t::eval_arg1_t;                                                     \
    using eval_arg2_t = typename base_t::eval_arg2_t;                                                     \
    using eval_t      = typename base_t::eval_t;                                                          \
                                                                                                          \
    name##FAExpression (const EArg1& left, const EArg2& right)                                            \
      : base_t(left,right)                                                                                \
      , m_fv_arg1 (scream::constants::fill_value<eval_left_t>)                                            \
      , m_fv_arg2 (scream::constants::fill_value<eval_right_t>)                                           \
    {}                                                                                                    \
                                                                                                          \
    template<typename... Args>                                                                            \
    KOKKOS_INLINE_FUNCTION                                                                                \
    eval_t eval(Args... args) const {                                                                     \
      if constexpr (not expr_l) {                                                                         \
        return eval_impl(this->m_left,this->m_right.eval(args...));                                       \
      } else if constexpr (not expr_r) {                                                                  \
        return eval_impl(this->m_left.eval(args...),this->m_right);                                       \
      } else {                                                                                            \
        return eval_impl(this->m_left.eval(args...),this->m_right.eval(args...));                         \
      }                                                                                                   \
    }                                                                                                     \
                                                                                                          \
  protected:                                                                                              \
                                                                                                          \
    KOKKOS_INLINE_FUNCTION                                                                                \
    eval_t eval_impl (const eval_left_t& l, const eval_right_t& r) const {                                \
      if (l==m_fv_arg1 or r==m_fv_arg2)                                                                   \
        return m_fv;                                                                                      \
      else                                                                                                \
        return base_t::eval_impl(l,r);                                                                    \
    }                                                                                                     \
                                                                                                          \
    eval_left_t  m_fv_l;                                                                                  \
    eval_right_t m_fv_r;                                                                                  \
    eval_t       m_fv;                                                                                    \
  };                                                                                                      \
                                                                                                          \
  template<typename EArg1, typename EArg2>                                                                \
  struct is_expr<name##FAExrepssion<EArg1,EArg2>> : std::true_type {};                                    \
  template<typename EArg1, typename EArg2>                                                                \
  struct eval_return<name##FAExpression<EArg1,EArg2>> : eval_return<name##Expression<EArg1,EArg2>> {};    \
                                                                                                          \
  /* function to create a fill-aware expression on the fly */                                             \
  template<typename EArg1, typename EArg2>                                                                \
  std::enable_if_t<is_expr_v<EArg1> or is_expr_v<EArg2>,name##FAExpression<EArg1,EArg2>>                  \
  fa_##impl (const EArg1& arg1, const EArg2& arg2)                                                        \
  {                                                                                                       \
    return name##FAExpression<EArg1,EArg2>(arg1,arg2);                                                    \
  }

// Create classes
//FILL_AWARE_BINARY_MATH_EXPR(max,Max)
//FILL_AWARE_BINARY_MATH_EXPR(min,Min)

template<typename EArg1, typename EArg2>                                                                \
class MaxFAExpression : public MaxExpression<EArg1,EArg2> {                                       \
public:                                                                                                 \
  using base_t = MaxExpression<EArg1,EArg2>;                                                         \
                                                                                                        \
  static constexpr bool expr_l = is_expr_v<EArg1>;                                                      \
  static constexpr bool expr_r = is_expr_v<EArg2>;                                                      \
                                                                                                        \
  using eval_arg1_t = typename base_t::eval_arg1_t;                                                     \
  using eval_arg2_t = typename base_t::eval_arg2_t;                                                     \
  using eval_t      = typename base_t::eval_t;                                                          \
                                                                                                        \
  MaxFAExpression (const EArg1& left, const EArg2& right)                                            \
    : base_t(left,right)                                                                                \
    , m_fv_arg1 (scream::constants::fill_value<eval_left_t>)                                            \
    , m_fv_arg2 (scream::constants::fill_value<eval_right_t>)                                           \
  {}                                                                                                    \
                                                                                                        \
  template<typename... Args>                                                                            \
  KOKKOS_INLINE_FUNCTION                                                                                \
  eval_t eval(Args... args) const {                                                                     \
    if constexpr (not expr_l) {                                                                         \
      return eval_impl(this->m_left,this->m_right.eval(args...));                                       \
    } else if constexpr (not expr_r) {                                                                  \
      return eval_impl(this->m_left.eval(args...),this->m_right);                                       \
    } else {                                                                                            \
      return eval_impl(this->m_left.eval(args...),this->m_right.eval(args...));                         \
    }                                                                                                   \
  }                                                                                                     \
                                                                                                        \
protected:                                                                                              \
                                                                                                        \
  KOKKOS_INLINE_FUNCTION                                                                                \
  eval_t eval_impl (const eval_left_t& l, const eval_right_t& r) const {                                \
    if (l==m_fv_arg1 or r==m_fv_arg2)                                                                   \
      return m_fv;                                                                                      \
    else                                                                                                \
      return base_t::eval_impl(l,r);                                                                    \
  }                                                                                                     \
                                                                                                        \
  eval_left_t  m_fv_l;                                                                                  \
  eval_right_t m_fv_r;                                                                                  \
  eval_t       m_fv;                                                                                    \
};                                                                                                      \
                                                                                                        \
template<typename EArg1, typename EArg2>                                                                \
struct is_expr<MaxFAExrepssion<EArg1,EArg2>> : std::true_type {};                                    \
template<typename EArg1, typename EArg2>                                                                \
struct eval_return<MaxFAExpression<EArg1,EArg2>> : eval_return<MaxExpression<EArg1,EArg2>> {};    \
                                                                                                        \
/* function to create a fill-aware expression on the fly */                                             \
template<typename EArg1, typename EArg2>                                                                \
std::enable_if_t<is_expr_v<EArg1> or is_expr_v<EArg2>,MaxFAExpression<EArg1,EArg2>>                  \
fa_max (const EArg1& arg1, const EArg2& arg2)                                                        \
{                                                                                                       \
  return MaxFAExpression<EArg1,EArg2>(arg1,arg2);                                                    \
}

#undef FILL_AWARE_BINARY_MATH_EXPR

} // namespace ekat

#endif // EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP
