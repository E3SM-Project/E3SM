#ifndef SCREAM_FIELD_IMPL_DETAILS_HPP
#define SCREAM_FIELD_IMPL_DETAILS_HPP

#include <ekat/kokkos/ekat_kokkos_types.hpp>

#include <vector>

namespace scream
{

namespace details {

template<CombineMode CM, bool use_fill, typename LhsView, typename RhsView, typename ST>
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
    if constexpr (use_fill)
      if constexpr (N==0)
        combine_and_fill<CM>(rhs(),lhs(),fill_val,alpha,beta);
      else
        combine_and_fill<CM>(rhs(i),lhs(i),fill_val,alpha,beta);
    else
      if constexpr (N==0)
        combine<CM>(rhs(),lhs(),alpha,beta);
      else
        combine<CM>(rhs(i),lhs(i),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    if constexpr (use_fill)
      combine_and_fill<CM>(rhs(i,j),lhs(i,j),fill_val,alpha,beta);
    else
      combine<CM>(rhs(i,j),lhs(i,j),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    if constexpr (use_fill)
      combine_and_fill<CM>(rhs(i,j,k),lhs(i,j,k),fill_val,alpha,beta);
    else
      combine<CM>(rhs(i,j,k),lhs(i,j,k),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    if constexpr (use_fill)
      combine_and_fill<CM>(rhs(i,j,k,l),lhs(i,j,k,l),fill_val,alpha,beta);
    else
      combine<CM>(rhs(i,j,k,l),lhs(i,j,k,l),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    if constexpr (use_fill)
      combine_and_fill<CM>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),fill_val,alpha,beta);
    else
      combine<CM>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    if constexpr (use_fill)
      combine_and_fill<CM>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),fill_val,alpha,beta);
    else
      combine<CM>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),alpha,beta);
  }

  ST alpha;
  ST beta;
  LhsView lhs;
  RhsView rhs;
  typename LhsView::traits::value_type fill_val;
};

template<CombineMode CM, bool use_fill, typename LhsView, typename RhsView, typename ST>
void
cvh (LhsView lhs, RhsView rhs,
     ST alpha, ST beta, typename LhsView::traits::value_type fill_val,
     const std::vector<int>& dims)
{
  CombineViewsHelper <CM, use_fill, LhsView, RhsView,  ST> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.fill_val = fill_val;
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

} // namespace scream

#endif // SCREAM_FIELD_IMPL_DETAILS_HPP
