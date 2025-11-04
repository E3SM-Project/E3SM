#ifndef EAMXX_MATH_FCNS_EXPRESSION_HPP
#define EAMXX_MATH_FCNS_EXPRESSION_HPP

#include "share/expressions/base.hpp"

namespace scream {

// ----------------- POW ------------------- //

template<typename EBase, typename EExp>
class PowExpression : public Expression<PowExpression<EBase,EExp>> {
public:
  static constexpr bool is_leaf = false;

  PowExpression (const EBase& base, const EExp& exp)
    : m_base(base)
    , m_exp(exp)
  {
    // Nothing to do here
  }

  int num_indices () const { return std::max(m_base.num_indices(),m_exp.num_indices()); }

  void set_eval_layout (const FieldLayout& fl) {
    m_base.set_eval_layout(fl);
    m_exp.set_eval_layout(fl);
  }

  KOKKOS_INLINE_FUNCTION
  Real eval () const {
    return Kokkos::pow(m_base.eval(),m_exp.eval());
  }
  KOKKOS_INLINE_FUNCTION
  void set_eval_data (const EvalData& data) const {
    m_base.set_eval_data(data);
    m_exp.set_eval_data(data);
  }
protected:

  EBase    m_base;
  EExp   m_exp;
};

template<typename EBase>
class PowExpression<EBase,Real> : public Expression<PowExpression<EBase,Real>> {
public:
  static constexpr bool is_leaf = false;

  PowExpression (const EBase& base, const Real exp)
    : m_base(base)
    , m_exp(exp)
  {
    // Nothing to do here
  }

  int num_indices () const { return m_base.num_indices(); }

  void set_eval_layout (const FieldLayout& fl) {
    m_base.set_eval_layout(fl);
  }

  KOKKOS_INLINE_FUNCTION
  Real eval () const {
    return Kokkos::pow(m_base.eval(),m_exp);
  }
  KOKKOS_INLINE_FUNCTION
  void set_eval_data (const EvalData& data) const {
    m_base.set_eval_data(data);
  }
protected:

  EBase    m_base;
  Real      m_exp;
};

template<typename EBase>
class PowExpression<EBase,int> : public Expression<PowExpression<EBase,int>> {
public:
  static constexpr bool is_leaf = false;

  PowExpression (const EBase& base, const int exp)
    : m_base(base)
    , m_exp(exp)
  {
    // Nothing to do here
  }

  int num_indices () const { return m_base.num_indices(); }

  void set_eval_layout (const FieldLayout& fl) {
    m_base.set_eval_layout(fl);
  }

  KOKKOS_INLINE_FUNCTION
  Real eval () const {
    return Kokkos::pow(m_base.eval(),m_exp);
  }
  KOKKOS_INLINE_FUNCTION
  void set_eval_data (const EvalData& data) const {
    m_base.set_eval_data(data);
  }
protected:

  EBase    m_base;
  int      m_exp;
};

template<typename EBase, typename EExp>
PowExpression<EBase,EExp>
pow (const Expression<EBase>& b, const Expression<EExp>& e)
{
  return PowExpression<EBase,EExp>(b.cast(),e.cast());
}

template<typename EBase>
PowExpression<EBase,int>
pow (const Expression<EBase>& b, const int e)
{
  return PowExpression<EBase,int>(b.cast(),e);
}

template<typename EBase>
PowExpression<EBase,Real>
pow (const Expression<EBase>& b, const Real e)
{
  return PowExpression<EBase,Real>(b.cast(),e);
}

// ----------------- Unary math fcns ------------------- //

#define UNARY_MATH_EXPRESSION(impl,name) \
  template<typename EArg>                                               \
  class name##Expression : public Expression<name##Expression<EArg>> {  \
  public:                                                               \
    name##Expression (const EArg& arg)                                  \
      : m_arg(arg)                                                      \
    {                                                                   \
      /* Nothing to do here */                                          \
    }                                                                   \
                                                                        \
    int num_indices () const { return m_arg.num_indices(); }            \
                                                                        \
    void set_eval_layout (const FieldLayout& fl) {                      \
      m_arg.set_eval_layout(fl);                                        \
    }                                                                   \
                                                                        \
    KOKKOS_INLINE_FUNCTION                                              \
    Real eval () const {                                                \
      return Kokkos::impl(m_arg.eval());                                \
    }                                                                   \
    KOKKOS_INLINE_FUNCTION                                              \
    void set_eval_data (const EvalData& data) const {                   \
      m_arg.set_eval_data(data);                                        \
    }                                                                   \
  protected:                                                            \
                                                                        \
    EArg    m_arg;                                                      \
  };                                                                    \
                                                                        \
  template<typename EArg>                                               \
  name##Expression<EArg>                                                \
  impl (const Expression<EArg>& arg)                                    \
  {                                                                     \
    return name##Expression<EArg>(arg.cast());                          \
  }

UNARY_MATH_EXPRESSION (sqrt,Sqrt)
UNARY_MATH_EXPRESSION (exp,Exp)
UNARY_MATH_EXPRESSION (log,Log)
UNARY_MATH_EXPRESSION (sin,Sin)
UNARY_MATH_EXPRESSION (cos,Cos)
#undef UNARY_MATH_EXPRESSION

} // namespace scream

#endif // EAMXX_MATH_FCNS_EXPRESSION_HPP
