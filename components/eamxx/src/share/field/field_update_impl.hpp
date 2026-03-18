#ifndef SCREAM_FIELD_UPDATE_IMPL_HPP
#define SCREAM_FIELD_UPDATE_IMPL_HPP

#include "share/field/field.hpp"

namespace scream
{

namespace details {

template<CombineMode CM, typename LhsView, typename RhsView, typename ST, typename MaskView, bool masked>
struct CombineViewsHelper {

  using exec_space = typename LhsView::traits::execution_space;

  static constexpr int N = LhsView::rank();
  static_assert( not masked or MaskView::rank()==N, "Mask view type has the wrong rank.\n" );

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
    if constexpr (masked) {
      if constexpr (N==0) {
        if (mask())
          combine<CM>(rhs(),lhs(),alpha,beta);
      } else {
        if (mask(i))
          combine<CM>(rhs(i),lhs(i),alpha,beta);
      }
    } else {
      if constexpr (N==0)
        combine<CM>(rhs(),lhs(),alpha,beta);
      else
        combine<CM>(rhs(i),lhs(i),alpha,beta);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    if constexpr (masked) {
      if (mask(i,j))
        combine<CM>(rhs(i,j),lhs(i,j),alpha,beta);
    } else {
      combine<CM>(rhs(i,j),lhs(i,j),alpha,beta);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    if constexpr (masked) {
      if (mask(i,j,k))
        combine<CM>(rhs(i,j,k),lhs(i,j,k),alpha,beta);
    } else {
      combine<CM>(rhs(i,j,k),lhs(i,j,k),alpha,beta);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    if constexpr (masked) {
      if (mask(i,j,k,l))
        combine<CM>(rhs(i,j,k,l),lhs(i,j,k,l),alpha,beta);
    } else {
      combine<CM>(rhs(i,j,k,l),lhs(i,j,k,l),alpha,beta);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    if constexpr (masked) {
      if (mask(i,j,k,l,m))
        combine<CM>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),alpha,beta);
    } else {
      combine<CM>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),alpha,beta);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    if constexpr (masked) {
      if (mask(i,j,k,l,m,n))
        combine<CM>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),alpha,beta);
    } else {
      combine<CM>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),alpha,beta);
    }
  }

  MaskView mask; 
  ST alpha;
  ST beta;
  LhsView lhs;
  RhsView rhs;
};

template<CombineMode CM, bool masked, typename LhsView, typename RhsView, typename MaskView, typename ST>
void
cvh (LhsView lhs, RhsView rhs,
     ST alpha, ST beta,
     const std::vector<int>& dims,
     MaskView mask)
{
  CombineViewsHelper <CM, LhsView, RhsView,  ST, MaskView, masked> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.mask = mask;
  helper.run(dims);
}

template<typename LhsView, typename MaskView, bool negate_mask>
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

  // In the following bool(mask(..)) returns true if mask!=0.
  // So doing bool(mask(..))!=negate_mask returns true if
  //  - mask!=0 and we don't need to negate the mask
  //  - mask==0 and we need to negate the mask

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if constexpr (N==0) {
      if (bool(mask())!=negate_mask)
        lhs() = value;
    } else {
      if (bool(mask(i))!=negate_mask)
        lhs(i) = value;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    if (bool(mask(i,j))!=negate_mask)
      lhs(i,j) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    if (bool(mask(i,j,k))!=negate_mask)
      lhs(i,j,k) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    if (bool(mask(i,j,k,l))!=negate_mask)
      lhs(i,j,k,l) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    if (bool(mask(i,j,k,l,m))!=negate_mask)
      lhs(i,j,k,l,m) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    if (bool(mask(i,j,k,l,m,n))!=negate_mask)
      lhs(i,j,k,l,m,n) = value;
  }

  ST value;
  LhsView lhs;
  MaskView mask;
};

template<bool negate_mask, typename LhsView, typename MaskView = LhsView>
void
svm (LhsView lhs,
     typename LhsView::traits::value_type value,
     const std::vector<int>& dims,
     MaskView mask = MaskView())
{
  SetValueMasked <LhsView, MaskView, negate_mask> helper;
  helper.lhs = lhs;
  helper.mask = mask;
  helper.value = value;

  EKAT_REQUIRE_MSG (mask.data()!=nullptr,
      "Error! Calling scream::details::svm with an invalid input mask view.\n");

  helper.run(dims);
}

} // namespace details

template<bool masked, CombineMode CM, typename ST, typename XST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta, const Field* mask)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where one of the two views is strided (or both)
  const auto x_contig = x.get_header().get_alloc_properties().contiguous();
  const auto y_contig = get_header().get_alloc_properties().contiguous();
  switch (layout.rank()) {
    case 0:
    {
      auto mv = masked ? mask->get_view<const int>() : view_dev_t<const int>{};
      if (x_contig and y_contig)
        details::cvh<CM,masked>(get_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims,mv);
      else if (x_contig)
        details::cvh<CM,masked>(get_strided_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims,mv);
      else if (y_contig)
        details::cvh<CM,masked>(get_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims,mv);
      else
        details::cvh<CM,masked>(get_strided_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims,mv);
      break;
    }
    case 1:
    {
      auto mv = masked ? mask->get_view<const int*>() : view_dev_t<const int*>{};
      if (x_contig and y_contig)
        details::cvh<CM,masked>(get_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims,mv);
      else if (x_contig)
        details::cvh<CM,masked>(get_strided_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims,mv);
      else if (y_contig)
        details::cvh<CM,masked>(get_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims,mv);
      else
        details::cvh<CM,masked>(get_strided_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims,mv);
      break;
    }
    case 2:
    {
      auto mv = masked ? mask->get_view<const int**>() : view_dev_t<const int**>{};
      if (x_contig and y_contig)
        details::cvh<CM,masked>(get_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims,mv);
      else if (x_contig)
        details::cvh<CM,masked>(get_strided_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims,mv);
      else if (y_contig)
        details::cvh<CM,masked>(get_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims,mv);
      else
        details::cvh<CM,masked>(get_strided_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims,mv);
      break;
    }
    case 3:
    {
      auto mv = masked ? mask->get_view<const int***>() : view_dev_t<const int***>{};
      if (x_contig and y_contig)
        details::cvh<CM,masked>(get_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims,mv);
      else if (x_contig)
        details::cvh<CM,masked>(get_strided_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims,mv);
      else if (y_contig)
        details::cvh<CM,masked>(get_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims,mv);
      else
        details::cvh<CM,masked>(get_strided_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims,mv);
      break;
    }
    case 4:
    {
      auto mv = masked ? mask->get_view<const int****>() : view_dev_t<const int****>{};
      if (x_contig and y_contig)
        details::cvh<CM,masked>(get_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims,mv);
      else if (x_contig)
        details::cvh<CM,masked>(get_strided_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims,mv);
      else if (y_contig)
        details::cvh<CM,masked>(get_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims,mv);
      else
        details::cvh<CM,masked>(get_strided_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims,mv);
      break;
    }
    case 5:
    {
      auto mv = masked ? mask->get_view<const int*****>() : view_dev_t<const int*****>{};
      if (x_contig and y_contig)
        details::cvh<CM,masked>(get_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims,mv);
      else if (x_contig)
        details::cvh<CM,masked>(get_strided_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims,mv);
      else if (y_contig)
        details::cvh<CM,masked>(get_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims,mv);
      else
        details::cvh<CM,masked>(get_strided_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims,mv);
      break;
    }
    case 6:
    {
      auto mv = masked ? mask->get_view<const int******>() : view_dev_t<const int******>{};
      if (x_contig and y_contig)
        details::cvh<CM,masked>(get_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims,mv);
      else if (x_contig)
        details::cvh<CM,masked>(get_strided_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims,mv);
      else if (y_contig)
        details::cvh<CM,masked>(get_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims,mv);
      else
        details::cvh<CM,masked>(get_strided_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims,mv);
      break;
    }
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_field.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n");
  }
  Kokkos::fence();
}

template<bool masked, bool negate_mask, typename ST>
void Field::deep_copy_impl (const ST value, const Field* mask)
{
  const auto& layout = get_header().get_identifier().get_layout();
  const auto  rank   = layout.rank();
  const auto& dims   = layout.dims();
  const auto contig = get_header().get_alloc_properties().contiguous();

  switch (rank) {
    case 0:
      if constexpr (masked) {
        if (contig)
          details::svm<negate_mask>(get_view<ST>(),value,dims,
                                    mask->get_view<const int>());
        else
          details::svm<negate_mask>(get_strided_view<ST>(),value,dims,
                                    mask->get_view<const int>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST>(),value);
      }
      break;
    case 1:
      if constexpr (masked) {
        if (contig)
          details::svm<negate_mask>(get_view<ST*>(),value,dims,
                                    mask->get_view<const int*>());
        else
          details::svm<negate_mask>(get_strided_view<ST*>(),value,dims,
                                    mask->get_view<const int*>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST*>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST*>(),value);
      }
      break;
    case 2:
      if constexpr (masked) {
        if (contig)
          details::svm<negate_mask>(get_view<ST**>(),value,dims,
                                    mask->get_view<const int**>());
        else
          details::svm<negate_mask>(get_strided_view<ST**>(),value,dims,
                                    mask->get_view<const int**>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST**>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST**>(),value);
      }
      break;
    case 3:
      if constexpr (masked) {
        if (contig)
          details::svm<negate_mask>(get_view<ST***>(),value,dims,
                                    mask->get_view<const int***>());
        else
          details::svm<negate_mask>(get_strided_view<ST***>(),value,dims,
                                    mask->get_view<const int***>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST***>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST***>(),value);
      }
      break;
    case 4:
      if constexpr (masked) {
        if (contig)
          details::svm<negate_mask>(get_view<ST****>(),value,dims,
                                    mask->get_view<const int****>());
        else
          details::svm<negate_mask>(get_strided_view<ST****>(),value,dims,
                                    mask->get_view<const int****>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST****>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST****>(),value);
      }
      break;
    case 5:
      if constexpr (masked) {
        if (contig)
          details::svm<negate_mask>(get_view<ST*****>(),value,dims,
                                    mask->get_view<const int*****>());
        else
          details::svm<negate_mask>(get_strided_view<ST*****>(),value,dims,
                                    mask->get_view<const int*****>());
      } else {
        if (contig)
          Kokkos::deep_copy(get_view<ST*****>(),value);
        else
          Kokkos::deep_copy(get_strided_view<ST*****>(),value);
      }
      break;
    case 6:
      if constexpr (masked) {
        if (contig)
          details::svm<negate_mask>(get_view<ST******>(),value,dims,
                                    mask->get_view<const int******>());
        else
          details::svm<negate_mask>(get_strided_view<ST******>(),value,dims,
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
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_IMPL_HPP
