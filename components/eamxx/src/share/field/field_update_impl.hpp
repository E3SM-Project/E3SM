#ifndef SCREAM_FIELD_UPDATE_IMPL_HPP
#define SCREAM_FIELD_UPDATE_IMPL_HPP

#include "share/field/field.hpp"

namespace scream
{

namespace details {

template<CombineMode CM, typename LhsView, typename RhsView, typename ST>
struct CombineViewsHelper {

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

  template<typename... Args>
  KOKKOS_INLINE_FUNCTION
  void operator() (Args... indices) const {
    auto& lhs_val = lhs.access(indices...);
    combine<CM>(rhs.access(indices...),lhs_val,alpha,beta);
    if constexpr (CM==CombineMode::Update)
      lhs_val += gamma;
  }

  ST alpha;
  ST beta;
  ST gamma;
  LhsView lhs;
  RhsView rhs;
};

template<CombineMode CM, typename LhsView, typename RhsView, typename ST>
void
cvh (LhsView lhs, RhsView rhs,
     ST alpha, ST beta, ST gamma,
     const std::vector<int>& dims)
{
  CombineViewsHelper <CM, LhsView, RhsView,  ST> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.gamma = gamma;
  helper.run(dims);
}

} // namespace details

template<CombineMode CM, typename ST, typename XST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta, const ST gamma)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where one of the two views is strided (or both)
  const auto x_lr_ok = x.get_header().get_alloc_properties().allows_layout_right();
  const auto y_lr_ok = get_header().get_alloc_properties().allows_layout_right();
  switch (layout.rank()) {
    case 0:
      if (x_lr_ok and y_lr_ok)
        details::cvh<CM>(get_view<ST>(),
                         x.get_view<const XST>(),
                         alpha,beta,gamma,dims);
      else if (x_lr_ok)
        details::cvh<CM>(get_strided_view<ST>(),
                         x.get_view<const XST>(),
                         alpha,beta,gamma,dims);
      else if (y_lr_ok)
        details::cvh<CM>(get_view<ST>(),
                         x.get_strided_view<const XST>(),
                         alpha,beta,gamma,dims);
      else
        details::cvh<CM>(get_strided_view<ST>(),
                         x.get_strided_view<const XST>(),
                         alpha,beta,gamma,dims);
      break;
    case 1:
      if (x_lr_ok and y_lr_ok)
        details::cvh<CM>(get_view<ST*>(),
                         x.get_view<const XST*>(),
                         alpha,beta,gamma,dims);
      else if (x_lr_ok)
        details::cvh<CM>(get_strided_view<ST*>(),
                         x.get_view<const XST*>(),
                         alpha,beta,gamma,dims);
      else if (y_lr_ok)
        details::cvh<CM>(get_view<ST*>(),
                         x.get_strided_view<const XST*>(),
                         alpha,beta,gamma,dims);
      else
        details::cvh<CM>(get_strided_view<ST*>(),
                         x.get_strided_view<const XST*>(),
                         alpha,beta,gamma,dims);
      break;
    case 2:
      if (x_lr_ok and y_lr_ok)
        details::cvh<CM>(get_view<ST**>(),
                         x.get_view<const XST**>(),
                         alpha,beta,gamma,dims);
      else if (x_lr_ok)
        details::cvh<CM>(get_strided_view<ST**>(),
                         x.get_view<const XST**>(),
                         alpha,beta,gamma,dims);
      else if (y_lr_ok)
        details::cvh<CM>(get_view<ST**>(),
                         x.get_strided_view<const XST**>(),
                         alpha,beta,gamma,dims);
      else
        details::cvh<CM>(get_strided_view<ST**>(),
                         x.get_strided_view<const XST**>(),
                         alpha,beta,gamma,dims);
      break;
    case 3:
      if (x_lr_ok and y_lr_ok)
        details::cvh<CM>(get_view<ST***>(),
                         x.get_view<const XST***>(),
                         alpha,beta,gamma,dims);
      else if (x_lr_ok)
        details::cvh<CM>(get_strided_view<ST***>(),
                         x.get_view<const XST***>(),
                         alpha,beta,gamma,dims);
      else if (y_lr_ok)
        details::cvh<CM>(get_view<ST***>(),
                         x.get_strided_view<const XST***>(),
                         alpha,beta,gamma,dims);
      else
        details::cvh<CM>(get_strided_view<ST***>(),
                         x.get_strided_view<const XST***>(),
                         alpha,beta,gamma,dims);
      break;
    case 4:
      if (x_lr_ok and y_lr_ok)
        details::cvh<CM>(get_view<ST****>(),
                         x.get_view<const XST****>(),
                         alpha,beta,gamma,dims);
      else if (x_lr_ok)
        details::cvh<CM>(get_strided_view<ST****>(),
                         x.get_view<const XST****>(),
                         alpha,beta,gamma,dims);
      else if (y_lr_ok)
        details::cvh<CM>(get_view<ST****>(),
                         x.get_strided_view<const XST****>(),
                         alpha,beta,gamma,dims);
      else
        details::cvh<CM>(get_strided_view<ST****>(),
                         x.get_strided_view<const XST****>(),
                         alpha,beta,gamma,dims);
      break;
    case 5:
      if (x_lr_ok and y_lr_ok)
        details::cvh<CM>(get_view<ST*****>(),
                         x.get_view<const XST*****>(),
                         alpha,beta,gamma,dims);
      else if (x_lr_ok)
        details::cvh<CM>(get_strided_view<ST*****>(),
                         x.get_view<const XST*****>(),
                         alpha,beta,gamma,dims);
      else if (y_lr_ok)
        details::cvh<CM>(get_view<ST*****>(),
                         x.get_strided_view<const XST*****>(),
                         alpha,beta,gamma,dims);
      else
        details::cvh<CM>(get_strided_view<ST*****>(),
                         x.get_strided_view<const XST*****>(),
                         alpha,beta,gamma,dims);
      break;
    case 6:
      if (x_lr_ok and y_lr_ok)
        details::cvh<CM>(get_view<ST******>(),
                         x.get_view<const XST******>(),
                         alpha,beta,gamma,dims);
      else if (x_lr_ok)
        details::cvh<CM>(get_strided_view<ST******>(),
                         x.get_view<const XST******>(),
                         alpha,beta,gamma,dims);
      else if (y_lr_ok)
        details::cvh<CM>(get_view<ST******>(),
                         x.get_strided_view<const XST******>(),
                         alpha,beta,gamma,dims);
      else
        details::cvh<CM>(get_strided_view<ST******>(),
                         x.get_strided_view<const XST******>(),
                         alpha,beta,gamma,dims);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_field.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n");
  }
  Kokkos::fence();
}

template<typename ST>
void Field::deep_copy_impl (const ST value)
{
  const auto rank = get_header().get_identifier().get_layout().rank();
  const auto lr_ok = get_header().get_alloc_properties().allows_layout_right();

  switch (rank) {
    case 0:
      if (lr_ok)
        Kokkos::deep_copy(get_view<ST>(),value);
      else
        Kokkos::deep_copy(get_strided_view<ST>(),value);
      break;
    case 1:
      if (lr_ok)
        Kokkos::deep_copy(get_view<ST*>(),value);
      else
        Kokkos::deep_copy(get_strided_view<ST*>(),value);
      break;
    case 2:
      if (lr_ok)
        Kokkos::deep_copy(get_view<ST**>(),value);
      else
        Kokkos::deep_copy(get_strided_view<ST**>(),value);
      break;
    case 3:
      if (lr_ok)
        Kokkos::deep_copy(get_view<ST***>(),value);
      else
        Kokkos::deep_copy(get_strided_view<ST***>(),value);
      break;
    case 4:
      if (lr_ok)
        Kokkos::deep_copy(get_view<ST****>(),value);
      else
        Kokkos::deep_copy(get_strided_view<ST****>(),value);
      break;
    case 5:
      if (lr_ok)
        Kokkos::deep_copy(get_view<ST*****>(),value);
      else
        Kokkos::deep_copy(get_strided_view<ST*****>(),value);
      break;
    case 6:
      if (lr_ok)
        Kokkos::deep_copy(get_view<ST******>(),value);
      else
        Kokkos::deep_copy(get_strided_view<ST******>(),value);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }
}


} // namespace scream

#endif // SCREAM_FIELD_UPDATE_IMPL_HPP
