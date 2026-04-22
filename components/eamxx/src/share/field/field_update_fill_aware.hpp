#ifndef SCREAM_FIELD_UPDATE_FILL_AWARE_HPP
#define SCREAM_FIELD_UPDATE_FILL_AWARE_HPP

#include "share/field/field.hpp"

namespace scream
{

namespace details {

template<CombineMode CM, typename LhsView, typename RhsView, typename ST>
struct CombineViewsFillAwareHelper {

  using exec_space = typename LhsView::traits::execution_space;

  static constexpr int N = LhsView::rank();

  template<int M>
  using MDRange = Kokkos::MDRangePolicy<
                    exec_space,
                    Kokkos::Rank<M,Kokkos::Iterate::Right,Kokkos::Iterate::Right>
                  >;

  void run (const std::vector<int>& dims) const {
    if constexpr (N==0) {
      Kokkos::RangePolicy<exec_space> policy(0,1);
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==1) {
      Kokkos::RangePolicy<exec_space> policy(0,dims[0]);
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==2) {
      MDRange<2> policy({0,0},{dims[0],dims[1]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==3) {
      MDRange<3> policy({0,0,0},{dims[0],dims[1],dims[2]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==4) {
      MDRange<4> policy({0,0,0,0},{dims[0],dims[1],dims[2],dims[3]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==5) {
      MDRange<5> policy({0,0,0,0,0},{dims[0],dims[1],dims[2],dims[3],dims[4]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==6) {
      MDRange<6> policy({0,0,0,0,0,0},{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]});
      Kokkos::parallel_for(policy,*this);
    } else {
      EKAT_ERROR_MSG ("Unsupported rank! Should be in [0,6].\n");
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if constexpr (N==0)
      combine_fill_aware<CM>(rhs(),lhs(),alpha,beta);
    else
      combine_fill_aware<CM>(rhs(i),lhs(i),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    combine_fill_aware<CM>(rhs(i,j),lhs(i,j),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    combine_fill_aware<CM>(rhs(i,j,k),lhs(i,j,k),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    combine_fill_aware<CM>(rhs(i,j,k,l),lhs(i,j,k,l),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    combine_fill_aware<CM>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    combine_fill_aware<CM>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),alpha,beta);
  }

  ST alpha;
  ST beta;
  LhsView lhs;
  RhsView rhs;
};

template<CombineMode CM, typename LhsView, typename RhsView, typename ST>
void
cvfah (LhsView lhs, RhsView rhs,
     ST alpha, ST beta,
     const std::vector<int>& dims)
{
  CombineViewsFillAwareHelper <CM, LhsView, RhsView,  ST> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.run(dims);
}

} // namespace details

template<CombineMode CM, typename ST, typename XST>
void Field::
update_fill_aware (const Field& x, const ST alpha, const ST beta)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where one of the two views is strided (or both)
  const auto x_lr_ok = x.get_header().get_alloc_properties().allows_layout_right();
  const auto y_lr_ok = get_header().get_alloc_properties().allows_layout_right();
  switch (layout.rank()) {
    case 0:
      if (x_lr_ok and y_lr_ok)
        details::cvfah<CM>(get_view<ST>(),
                           x.get_view<const XST>(),
                           alpha,beta,dims);
      else if (x_lr_ok)
        details::cvfah<CM>(get_strided_view<ST>(),
                           x.get_view<const XST>(),
                           alpha,beta,dims);
      else if (y_lr_ok)
        details::cvfah<CM>(get_view<ST>(),
                           x.get_strided_view<const XST>(),
                           alpha,beta,dims);
      else
        details::cvfah<CM>(get_strided_view<ST>(),
                           x.get_strided_view<const XST>(),
                           alpha,beta,dims);
      break;
    case 1:
      if (x_lr_ok and y_lr_ok)
        details::cvfah<CM>(get_view<ST*>(),
                           x.get_view<const XST*>(),
                           alpha,beta,dims);
      else if (x_lr_ok)
        details::cvfah<CM>(get_strided_view<ST*>(),
                           x.get_view<const XST*>(),
                           alpha,beta,dims);
      else if (y_lr_ok)
        details::cvfah<CM>(get_view<ST*>(),
                           x.get_strided_view<const XST*>(),
                           alpha,beta,dims);
      else
        details::cvfah<CM>(get_strided_view<ST*>(),
                           x.get_strided_view<const XST*>(),
                           alpha,beta,dims);
      break;
    case 2:
      if (x_lr_ok and y_lr_ok)
        details::cvfah<CM>(get_view<ST**>(),
                           x.get_view<const XST**>(),
                           alpha,beta,dims);
      else if (x_lr_ok)
        details::cvfah<CM>(get_strided_view<ST**>(),
                           x.get_view<const XST**>(),
                           alpha,beta,dims);
      else if (y_lr_ok)
        details::cvfah<CM>(get_view<ST**>(),
                           x.get_strided_view<const XST**>(),
                           alpha,beta,dims);
      else
        details::cvfah<CM>(get_strided_view<ST**>(),
                           x.get_strided_view<const XST**>(),
                           alpha,beta,dims);
      break;
    case 3:
      if (x_lr_ok and y_lr_ok)
        details::cvfah<CM>(get_view<ST***>(),
                           x.get_view<const XST***>(),
                           alpha,beta,dims);
      else if (x_lr_ok)
        details::cvfah<CM>(get_strided_view<ST***>(),
                           x.get_view<const XST***>(),
                           alpha,beta,dims);
      else if (y_lr_ok)
        details::cvfah<CM>(get_view<ST***>(),
                           x.get_strided_view<const XST***>(),
                           alpha,beta,dims);
      else
        details::cvfah<CM>(get_strided_view<ST***>(),
                           x.get_strided_view<const XST***>(),
                           alpha,beta,dims);
      break;
    case 4:
      if (x_lr_ok and y_lr_ok)
        details::cvfah<CM>(get_view<ST****>(),
                           x.get_view<const XST****>(),
                           alpha,beta,dims);
      else if (x_lr_ok)
        details::cvfah<CM>(get_strided_view<ST****>(),
                           x.get_view<const XST****>(),
                           alpha,beta,dims);
      else if (y_lr_ok)
        details::cvfah<CM>(get_view<ST****>(),
                           x.get_strided_view<const XST****>(),
                           alpha,beta,dims);
      else
        details::cvfah<CM>(get_strided_view<ST****>(),
                           x.get_strided_view<const XST****>(),
                           alpha,beta,dims);
      break;
    case 5:
      if (x_lr_ok and y_lr_ok)
        details::cvfah<CM>(get_view<ST*****>(),
                           x.get_view<const XST*****>(),
                           alpha,beta,dims);
      else if (x_lr_ok)
        details::cvfah<CM>(get_strided_view<ST*****>(),
                           x.get_view<const XST*****>(),
                           alpha,beta,dims);
      else if (y_lr_ok)
        details::cvfah<CM>(get_view<ST*****>(),
                           x.get_strided_view<const XST*****>(),
                           alpha,beta,dims);
      else
        details::cvfah<CM>(get_strided_view<ST*****>(),
                           x.get_strided_view<const XST*****>(),
                           alpha,beta,dims);
      break;
    case 6:
      if (x_lr_ok and y_lr_ok)
        details::cvfah<CM>(get_view<ST******>(),
                           x.get_view<const XST******>(),
                           alpha,beta,dims);
      else if (x_lr_ok)
        details::cvfah<CM>(get_strided_view<ST******>(),
                           x.get_view<const XST******>(),
                           alpha,beta,dims);
      else if (y_lr_ok)
        details::cvfah<CM>(get_view<ST******>(),
                           x.get_strided_view<const XST******>(),
                           alpha,beta,dims);
      else
        details::cvfah<CM>(get_strided_view<ST******>(),
                           x.get_strided_view<const XST******>(),
                           alpha,beta,dims);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_field.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n");
  }
  Kokkos::fence();
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_FILL_AWARE_HPP
