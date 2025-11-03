#ifndef EAMXX_FIELD_EXPRESSION_HPP
#define EAMXX_FIELD_EXPRESSION_HPP

#include "share/expressions/base.hpp"

#include "share/field/field.hpp"

namespace scream {

// TODO: support 4+ dim. Also 0d?
template<typename DT>
class FieldExpressionBase : public Expression<FieldExpressionBase<DT>> {
public:
  static constexpr bool is_leaf = true;

  FieldExpressionBase (const Field& f)
  {
    EKAT_REQUIRE_MSG (get_data_type<DT>()==f.data_type(),
        "Error! Wrong template arg for this field.\n");
    switch (f.rank()) {
      case 1: m_view_1d = f.get_view<const DT*>();   break;
      case 2: m_view_2d = f.get_view<const DT**>();  break;
      case 3: m_view_3d = f.get_view<const DT***>(); break;
      default:
        EKAT_ERROR_MSG ("Unsupported rank (" << f.rank() << ").\n");
    }

    m_rank = f.rank();
  }

  int num_indices () const { return m_rank; }

  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval(int i) const {
    return static_cast<T>(m_view_1d(i));
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval(int i, int j) const {
    return static_cast<T>(m_view_2d(i,j));
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval(int i,int j,int k) const {
    return static_cast<T>(m_view_3d(i,j,k));
  }

  // void set_eval_layout (const FieldLayout& fl) {}
protected:

  using KT = typename Field::kt_dev;

  typename KT::view_1d<const DT> m_view_1d;
  typename KT::view_2d<const DT> m_view_2d;
  typename KT::view_3d<const DT> m_view_3d;

  int m_rank;
};

using RealFieldExpression = FieldExpressionBase<Real>;
using IntFieldExpression  = FieldExpressionBase<int>;

} // namespace scream

#endif // EAMXX_FIELD_EXPRESSION_HPP
