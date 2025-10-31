#ifndef EAMXX_FIELD_EXPRESSION_HPP
#define EAMXX_FIELD_EXPRESSION_HPP

#include "share/expressions/base.hpp"

#include "share/field/field.hpp"

namespace scream {

// TODO: support 4+ dim. Also 0d?
template<typename DT>
class FieldExpressionBase : public Expression<FieldExpressionBase<DT>,DT> {
public:
  static constexpr bool is_leaf = true;

  using ret_t = DT;

  FieldExpressionBase (const Field& f)
  {
    EKAT_REQUIRE_MSG (get_data_type<DT>()==f.data_type(),
        "Error! Wrong template arg for this field.\n");
    switch (f.rank()) {
      case 1: m_view_1d = f.get_view<const ret_t*>();   break;
      case 2: m_view_2d = f.get_view<const ret_t**>();  break;
      case 3: m_view_3d = f.get_view<const ret_t***>(); break;
      default:
        EKAT_ERROR_MSG ("Unsupported rank (" << f.rank() << ").\n");
    }
  }

  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i) const {
    return m_view_1d(i);
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i, int j) const {
    return m_view_2d(i,j);
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i, int j, int k) const {
    return m_view_3d(i,j,k);
  }

protected:

  using KT = typename Field::kt_dev;

  typename KT::view_1d<const ret_t> m_view_1d;
  typename KT::view_2d<const ret_t> m_view_2d;
  typename KT::view_3d<const ret_t> m_view_3d;
};

using RealFieldExpression = FieldExpressionBase<Real>;
using IntFieldExpression  = FieldExpressionBase<int>;

} // namespace scream

#endif // EAMXX_FIELD_EXPRESSION_HPP
