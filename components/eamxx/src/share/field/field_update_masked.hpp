#ifndef SCREAM_FIELD_UPDATE_MASKED_HPP
#define SCREAM_FIELD_UPDATE_MASKED_HPP

#include "share/field/field.hpp"
#include "share/field/field_dispatch.hpp"

namespace scream
{

namespace details {

template<CombineMode CM, typename LhsView, typename RhsView, typename ST, typename MaskView>
struct CombineViewsMaskedHelper {

  using exec_space = typename LhsView::traits::execution_space;

  static constexpr int N = LhsView::rank();
  static_assert( MaskView::rank()==N, "Mask view type has the wrong rank.\n" );

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
    if (mask.access(indices...)) {
      auto& lhs_val = lhs.access(indices...);
      combine<CM>(rhs.access(indices...),lhs_val,alpha,beta);
      if constexpr (CM==CombineMode::Update)
        lhs_val += gamma;
    }
  }

  MaskView mask; 
  ST alpha;
  ST beta;
  ST gamma;
  LhsView lhs;
  RhsView rhs;
};

template<CombineMode CM, typename LhsView, typename RhsView, typename MaskView, typename ST>
void
cvmh (LhsView lhs, RhsView rhs,
      ST alpha, ST beta, ST gamma,
      const std::vector<int>& dims,
      MaskView mask)
{
  CombineViewsMaskedHelper <CM, LhsView, RhsView,  ST, MaskView> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.gamma = gamma;
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
      EKAT_ERROR_MSG ("Unsupported rank! Should be in [0,6].\n");
    }
  }

  template<typename... Args>
  KOKKOS_INLINE_FUNCTION
  void operator() (Args... indices) const {
    // Execute block if mask=0 and negate_mask=true or mask=1 and negate_mask=false
    if ( (mask.access(indices...)!=0) != negate_mask )
      lhs.access(indices...) = value;
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
     MaskView mask)
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

template<CombineMode CM, typename ST, typename XST>
void Field::
update_masked (const Field& x, const ST alpha, const ST beta, const ST gamma, const Field& mask) const
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where any of the views is strided
  details::dispatch_views(layout.rank(),
    [&](const auto& y_view, const auto& x_view, const auto& m_view) {
      details::cvmh<CM>(y_view,x_view,alpha,beta,gamma,dims,m_view);
    },
    details::as<ST>(*this), details::as<const XST>(x), details::as<const int>(mask));
  Kokkos::fence();
}

template<bool negate_mask, typename ST>
void Field::deep_copy_masked (const ST value, const Field& mask) const
{
  const auto& layout = get_header().get_identifier().get_layout();
  const auto& dims   = layout.dims();

  details::dispatch_views(layout.rank(),
    [&](const auto& y_view, const auto& m_view) {
      details::svm<negate_mask>(y_view,value,dims,m_view);
    },
    details::as<ST>(*this), details::as<const int>(mask));
  Kokkos::fence();
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_MASKED_HPP
