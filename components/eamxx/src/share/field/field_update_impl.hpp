#ifndef SCREAM_FIELD_UPDATE_IMPL_HPP
#define SCREAM_FIELD_UPDATE_IMPL_HPP

#include "share/field/field.hpp"

namespace scream
{

namespace details {

template<CombineMode CM, bool FillAware, typename LhsView, typename RhsView, typename ST>
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
      combine<CM,FillAware>(rhs(),lhs(),alpha,beta);
    else
      combine<CM,FillAware>(rhs(i),lhs(i),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    combine<CM,FillAware>(rhs(i,j),lhs(i,j),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    combine<CM,FillAware>(rhs(i,j,k),lhs(i,j,k),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    combine<CM,FillAware>(rhs(i,j,k,l),lhs(i,j,k,l),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    combine<CM,FillAware>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    combine<CM,FillAware>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),alpha,beta);
  }

  ST alpha;
  ST beta;
  LhsView lhs;
  RhsView rhs;
};

template<CombineMode CM, bool FillAware, typename LhsView, typename RhsView, typename ST>
void
cvh (LhsView lhs, RhsView rhs,
     ST alpha, ST beta,
     const std::vector<int>& dims)
{
  CombineViewsHelper <CM, FillAware, LhsView, RhsView,  ST> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.run(dims);
}

template<typename LhsView, typename MaskView, bool use_mask>
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
      lhs_ref = not use_mask or mask() ? value : lhs_ref;
    } else {
      auto& lhs_ref = lhs(i);
      lhs_ref = not use_mask or mask(i) ? value : lhs_ref;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    auto& lhs_ref = lhs(i,j);
    lhs_ref = not use_mask or mask(i,j) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    auto& lhs_ref = lhs(i,j,k);
    lhs_ref = not use_mask or mask(i,j,k) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    auto& lhs_ref = lhs(i,j,k,l);
    lhs_ref = not use_mask or mask(i,j,k,l) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    auto& lhs_ref = lhs(i,j,k,l,m);
    lhs_ref = not use_mask or mask(i,j,k,l,m) ? value : lhs_ref;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    auto& lhs_ref = lhs(i,j,k,l,m,n);
    lhs_ref = not use_mask or mask(i,j,k,l,m,n) ? value : lhs_ref;
  }

  ST value;
  LhsView lhs;
  MaskView mask;
};

template<bool use_mask, typename LhsView, typename MaskView = LhsView>
void
svm (LhsView lhs,
     typename LhsView::traits::value_type value,
     const std::vector<int>& dims,
     MaskView mask = MaskView())
{
  SetValueMasked <LhsView, MaskView, use_mask> helper;
  helper.lhs = lhs;
  helper.mask = mask;
  helper.value = value;

  EKAT_REQUIRE_MSG (not use_mask or mask.data()!=nullptr,
      "Error! Calling scream::details::svm with use_mask=true, but input mask view is invalid.\n");

  helper.run(dims);
}

} // namespace details

template<CombineMode CM, bool FillAware, typename ST, typename XST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where one of the two views is strided (or both)
  const auto x_contig = x.get_header().get_alloc_properties().contiguous();
  const auto y_contig = get_header().get_alloc_properties().contiguous();
  switch (layout.rank()) {
    case 0:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims);
      break;
    case 1:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims);
      break;
    case 2:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims);
      break;
    case 3:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims);
      break;
    case 4:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims);
      break;
    case 5:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims);
      break;
    case 6:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST******>(),
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

template<bool use_mask, typename ST>
void Field::deep_copy_impl (const ST value, const Field& mask)
{
  const auto& layout = get_header().get_identifier().get_layout();
  const auto  rank   = layout.rank();
  const auto& dims   = layout.dims();
  const auto contig = get_header().get_alloc_properties().contiguous();

  switch (rank) {
    case 0:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST>(),value,dims,
                                 mask.get_view<const int>());
        else
          details::svm<use_mask>(get_strided_view<ST>(),value,dims,
                                 mask.get_view<const int>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST>(),value,dims);
      }
      break;
    case 1:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST*>(),value,dims,
                                 mask.get_view<const int*>());
        else
          details::svm<use_mask>(get_strided_view<ST*>(),value,dims,
                                 mask.get_view<const int*>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST*>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST*>(),value,dims);
      }
      break;
    case 2:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST**>(),value,dims,
                                 mask.get_view<const int**>());
        else
          details::svm<use_mask>(get_strided_view<ST**>(),value,dims,
                                 mask.get_view<const int**>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST**>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST**>(),value,dims);
      }
      break;
    case 3:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST***>(),value,dims,
                                 mask.get_view<const int***>());
        else
          details::svm<use_mask>(get_strided_view<ST***>(),value,dims,
                                 mask.get_view<const int***>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST***>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST***>(),value,dims);
      }
      break;
    case 4:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST****>(),value,dims,
                                 mask.get_view<const int****>());
        else
          details::svm<use_mask>(get_strided_view<ST****>(),value,dims,
                                 mask.get_view<const int****>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST****>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST****>(),value,dims);
      }
      break;
    case 5:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST*****>(),value,dims,
                                 mask.get_view<const int*****>());
        else
          details::svm<use_mask>(get_strided_view<ST*****>(),value,dims,
                                 mask.get_view<const int*****>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST*****>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST*****>(),value,dims);
      }
      break;
    case 6:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST******>(),value,dims,
                                 mask.get_view<const int******>());
        else
          details::svm<use_mask>(get_strided_view<ST******>(),value,dims,
                                 mask.get_view<const int******>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST******>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST******>(),value,dims);
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_IMPL_HPP
