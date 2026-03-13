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
      EKAT_ERROR_MSG ("Unsupported rank! Should be in [2,6].\n");
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if constexpr (N==0)
      combine<CM>(rhs(),lhs(),alpha,beta,fvh);
    else
      combine<CM>(rhs(i),lhs(i),alpha,beta,fvh);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    combine<CM>(rhs(i,j),lhs(i,j),alpha,beta,fvh);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    combine<CM>(rhs(i,j,k),lhs(i,j,k),alpha,beta,fvh);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    combine<CM>(rhs(i,j,k,l),lhs(i,j,k,l),alpha,beta,fvh);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    combine<CM>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),alpha,beta,fvh);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    combine<CM>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),alpha,beta,fvh);
  }

  FillValueHandling fvh;
  ST alpha;
  ST beta;
  LhsView lhs;
  RhsView rhs;
};

template<CombineMode CM, typename LhsView, typename RhsView, typename ST>
void
cvh (LhsView lhs, RhsView rhs,
     ST alpha, ST beta,
     const std::vector<int>& dims,
     FillValueHandling fvh = FillValueHandling::None)
{
  CombineViewsHelper <CM, LhsView, RhsView,  ST> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.fvh = fvh;
  helper.run(dims);
}

template<typename LhsView, typename MaskView>
struct SetValueMasked
{
  using exec_space = typename LhsView::traits::execution_space;

  using ST = typename LhsView::traits::value_type;

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
      EKAT_ERROR_MSG ("Unsupported rank! Should be in [2,6].\n");
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if constexpr (N==0) {
      auto& lhs_ref = lhs();
      lhs_ref = mask() ? value : lhs_ref;
    } else {
      auto& lhs_ref = lhs(i);
      lhs_ref = mask(i) ? value : lhs_ref;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    auto& lhs_ref = lhs(i,j);
    lhs_ref = mask(i,j) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    auto& lhs_ref = lhs(i,j,k);
    lhs_ref = mask(i,j,k) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    auto& lhs_ref = lhs(i,j,k,l);
    lhs_ref = mask(i,j,k,l) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    auto& lhs_ref = lhs(i,j,k,l,m);
    lhs_ref = mask(i,j,k,l,m) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    auto& lhs_ref = lhs(i,j,k,l,m,n);
    lhs_ref = mask(i,j,k,l,m,n) ? value : lhs_ref;
  }

  ST value;
  LhsView lhs;
  MaskView mask;
};

template<typename LhsView, typename MaskView = LhsView>
void
svm (LhsView lhs,
     typename LhsView::traits::value_type value,
     const std::vector<int>& dims,
     MaskView mask = MaskView())
{
  SetValueMasked <LhsView, MaskView> helper;
  helper.lhs = lhs;
  helper.mask = mask;
  helper.value = value;

  EKAT_REQUIRE_MSG (mask.data()!=nullptr,
      "Error! Calling scream::details::svm with an invalid input mask view.\n");

  helper.run(dims);
}

} // namespace details

template<CombineMode CM, typename ST, typename XST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Check how we need to handle fill_value entries.
  // If we know that RHS won't contain fill_value entries AND this field was
  // set to ignore fv entries in the RHS, we can override fvh with None
  // NOTE: Absorbing needs to remain untouched, since LHS may have fill_values.
  // NOTE: Replace ignores fill value handling, so y=x REGARDLESS of x's content
  auto fvh = get_header().fill_value_handling();
  auto x_may_have_fv = x.get_header().may_be_filled();
  if (CM==CombineMode::Replace or (not x.get_header().may_be_filled() and fvh==IgnoreRhs)) {
    fvh = None;
  }

  // If the LHS is "absorbing" FV entries AND the rhs MAY have FV entries,
  // we MUST flag the LHS as potentially filled
  if (fvh==Absorbing and x_may_have_fv) {
    get_header().set_may_be_filled(true);
  }

  // Must handle the case where one of the two views is strided (or both)
  const auto x_contig = x.get_header().get_alloc_properties().contiguous();
  const auto y_contig = get_header().get_alloc_properties().contiguous();
  switch (layout.rank()) {
    case 0:
      if (x_contig and y_contig)
        details::cvh<CM>(get_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims,fvh);
      else if (x_contig)
        details::cvh<CM>(get_strided_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims,fvh);
      else if (y_contig)
        details::cvh<CM>(get_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims,fvh);
      else
        details::cvh<CM>(get_strided_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims,fvh);
      break;
    case 1:
      if (x_contig and y_contig)
        details::cvh<CM>(get_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims,fvh);
      else if (x_contig)
        details::cvh<CM>(get_strided_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims,fvh);
      else if (y_contig)
        details::cvh<CM>(get_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims,fvh);
      else
        details::cvh<CM>(get_strided_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims,fvh);
      break;
    case 2:
      if (x_contig and y_contig)
        details::cvh<CM>(get_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims,fvh);
      else if (x_contig)
        details::cvh<CM>(get_strided_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims,fvh);
      else if (y_contig)
        details::cvh<CM>(get_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims,fvh);
      else
        details::cvh<CM>(get_strided_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims,fvh);
      break;
    case 3:
      if (x_contig and y_contig)
        details::cvh<CM>(get_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims,fvh);
      else if (x_contig)
        details::cvh<CM>(get_strided_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims,fvh);
      else if (y_contig)
        details::cvh<CM>(get_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims,fvh);
      else
        details::cvh<CM>(get_strided_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims,fvh);
      break;
    case 4:
      if (x_contig and y_contig)
        details::cvh<CM>(get_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims,fvh);
      else if (x_contig)
        details::cvh<CM>(get_strided_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims,fvh);
      else if (y_contig)
        details::cvh<CM>(get_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims,fvh);
      else
        details::cvh<CM>(get_strided_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims,fvh);
      break;
    case 5:
      if (x_contig and y_contig)
        details::cvh<CM>(get_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims,fvh);
      else if (x_contig)
        details::cvh<CM>(get_strided_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims,fvh);
      else if (y_contig)
        details::cvh<CM>(get_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims,fvh);
      else
        details::cvh<CM>(get_strided_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims,fvh);
      break;
    case 6:
      if (x_contig and y_contig)
        details::cvh<CM>(get_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims,fvh);
      else if (x_contig)
        details::cvh<CM>(get_strided_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims,fvh);
      else if (y_contig)
        details::cvh<CM>(get_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims,fvh);
      else
        details::cvh<CM>(get_strided_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims,fvh);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_field.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n");
  }
  Kokkos::fence();
}

template<typename ST>
void Field::deep_copy_impl (const ST value, const Field* mask)
{
  const auto& layout = get_header().get_identifier().get_layout();
  const auto  rank   = layout.rank();
  const auto& dims   = layout.dims();
  const auto contig = get_header().get_alloc_properties().contiguous();

  switch (rank) {
    case 0:
      if (mask) {
        if (contig)
          details::svm(get_view<ST>(),value,dims,
                       mask->get_view<const int>());
        else
          details::svm(get_strided_view<ST>(),value,dims,
                       mask->get_view<const int>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST>(),value);
      }
      break;
    case 1:
      if (mask) {
        if (contig)
          details::svm(get_view<ST*>(),value,dims,
                       mask->get_view<const int*>());
        else
          details::svm(get_strided_view<ST*>(),value,dims,
                       mask->get_view<const int*>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST*>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST*>(),value);
      }
      break;
    case 2:
      if (mask) {
        if (contig)
          details::svm(get_view<ST**>(),value,dims,
                       mask->get_view<const int**>());
        else
          details::svm(get_strided_view<ST**>(),value,dims,
                       mask->get_view<const int**>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST**>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST**>(),value);
      }
      break;
    case 3:
      if (mask) {
        if (contig)
          details::svm(get_view<ST***>(),value,dims,
                       mask->get_view<const int***>());
        else
          details::svm(get_strided_view<ST***>(),value,dims,
                       mask->get_view<const int***>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST***>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST***>(),value);
      }
      break;
    case 4:
      if (mask) {
        if (contig)
          details::svm(get_view<ST****>(),value,dims,
                       mask->get_view<const int****>());
        else
          details::svm(get_strided_view<ST****>(),value,dims,
                       mask->get_view<const int****>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST****>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST****>(),value);
      }
      break;
    case 5:
      if (mask) {
        if (contig)
          details::svm(get_view<ST*****>(),value,dims,
                       mask->get_view<const int*****>());
        else
          details::svm(get_strided_view<ST*****>(),value,dims,
                       mask->get_view<const int*****>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST*****>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST*****>(),value);
      }
      break;
    case 6:
      if (mask) {
        if (contig)
          details::svm(get_view<ST******>(),value,dims,
                       mask->get_view<const int******>());
        else
          details::svm(get_strided_view<ST******>(),value,dims,
                       mask->get_view<const int******>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST******>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST******>(),value);
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }

  // If we're setting the field to fv, then it DEFINITELY could contain fv entries
  if (value==constants::fill_value<ST>) {
    get_header().set_may_be_filled(true);
  }
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_IMPL_HPP
