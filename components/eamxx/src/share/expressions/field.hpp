#ifndef EAMXX_FIELD_EXPRESSION_HPP
#define EAMXX_FIELD_EXPRESSION_HPP

#include "share/expressions/base.hpp"

#include "share/field/field.hpp"

#include <ekat_view_broadcast.hpp>

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

    const auto& fl = f.get_header().get_identifier().get_layout();
    std::copy_n(fl.tags().begin(),fl.rank(),m_field_tags);
    m_field_rank = f.rank();
  }

  // // FieldExpressionBase (const FieldExpressionBase& src) {
  // //   *this = src;
  // // }
  // FieldExpressionBase& operator= (const FieldExpressionBase&) = default;
  // FieldExpressionBase& operator= (const FieldExpressionBase& src) {
  //   m_view_1d = src.m_view_1d;
  //   m_view_2d = src.m_view_2d;
  //   m_view_3d = src.m_view_3d;

  //   m_eval_view_1d = src.m_eval_view_1d;
  //   m_eval_view_2d = src.m_eval_view_2d;
  //   m_eval_view_3d = src.m_eval_view_3d;

  //   m_data = src.m_data;

  //   // Cannot store FieldLayout, as it is not device friendly and would cause many warnings
  //   for (int i=0; i<src.m_field_rank; ++i)
  //     m_field_tags[i] = src.m_field_tags[i];

  //   m_field_rank = src.m_field_rank;
  //   m_eval_rank  = src.m_eval_rank;

  //   return *this;
  // }

  int num_indices () const { return m_field_rank; }

  KOKKOS_INLINE_FUNCTION
  Real eval() const {
    if (m_eval_rank==1)
      return m_eval_view_1d(m_data.i);
    else if (m_eval_rank==2)
      return m_eval_view_2d(m_data.i,m_data.j);
    else
      return m_eval_view_3d(m_data.i,m_data.j,m_data.k);
  }
  // KOKKOS_INLINE_FUNCTION
  // Real eval(int i) const {
  //   return static_cast<Real>(m_view_1d(i));
  // }
  // KOKKOS_INLINE_FUNCTION
  // Real eval(int i, int j) const {
  //   return static_cast<Real>(m_view_2d(i,j));
  // }
  // KOKKOS_INLINE_FUNCTION
  // Real eval(int i,int j,int k) const {
  //   return static_cast<Real>(m_view_3d(i,j,k));
  // }

  void set_eval_layout (const FieldLayout& fl) {
    m_eval_rank = fl.rank();
    EKAT_REQUIRE_MSG (m_eval_rank>=m_field_rank,
      "[FieldExpressionBase] Error! Evaluation rank must be higher than the field rank.\n"
      " - field rank: " << m_field_rank + "\n"
      " - eval rank : " << m_eval_rank + "\n");

    auto eval_dims = fl.dims();
    for (int i=0; i<m_field_rank; ++i) {
      int idx = fl.dim_idx(m_field_tags[i]);
      eval_dims[idx] = -1;
    }
    switch (m_eval_rank) {
      case 1:
        m_eval_view_1d.setup(m_view_1d,eval_dims);
        break;
      case 2:
        if (m_field_rank==1)
          m_eval_view_1d.setup(m_view_1d,eval_dims);
        else
          m_eval_view_1d.setup(m_view_2d,eval_dims);
        break;
      case 3:
        if (m_field_rank==1)
          m_eval_view_1d.setup(m_view_1d,eval_dims);
        else if (m_field_rank==1)
          m_eval_view_1d.setup(m_view_2d,eval_dims);
        else
          m_eval_view_1d.setup(m_view_3d,eval_dims);
        break;
      default:
        EKAT_ERROR_MSG (
          "[FieldExpressionBase] Error! Unsupported field evaluation rank.\n"
          " - eval rank : " << m_eval_rank << "\n");
    }
  }
  KOKKOS_INLINE_FUNCTION
  void set_eval_data (const EvalData& data) const { m_data = data; }
protected:

  using KT = typename Field::kt_dev;
  template<typename T, int N>
  using view_nd = typename KT::view_ND<T,N>;

  view_nd<const DT,1> m_view_1d;
  view_nd<const DT,2> m_view_2d;
  view_nd<const DT,3> m_view_3d;

  ekat::ViewBroadcast<view_nd<const DT,1>> m_eval_view_1d;
  ekat::ViewBroadcast<view_nd<const DT,2>> m_eval_view_2d;
  ekat::ViewBroadcast<view_nd<const DT,3>> m_eval_view_3d;

  mutable EvalData m_data;

  // Cannot store FieldLayout, as it is not device friendly and would cause many warnings
  FieldTag m_field_tags[8];

  int m_field_rank;
  int m_eval_rank;
};

using RealFieldExpression = FieldExpressionBase<Real>;
using IntFieldExpression  = FieldExpressionBase<int>;

} // namespace scream

#endif // EAMXX_FIELD_EXPRESSION_HPP
