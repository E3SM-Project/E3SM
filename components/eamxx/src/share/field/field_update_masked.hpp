#ifndef SCREAM_FIELD_UPDATE_MASKED_HPP
#define SCREAM_FIELD_UPDATE_MASKED_HPP

#include "share/field/field.hpp"

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

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if constexpr (N==0) {
      if (mask())
        combine<CM>(rhs(),lhs(),alpha,beta);
    } else {
      if (mask(i))
        combine<CM>(rhs(i),lhs(i),alpha,beta);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    if (mask(i,j))
      combine<CM>(rhs(i,j),lhs(i,j),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    if (mask(i,j,k))
      combine<CM>(rhs(i,j,k),lhs(i,j,k),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    if (mask(i,j,k,l))
      combine<CM>(rhs(i,j,k,l),lhs(i,j,k,l),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    if (mask(i,j,k,l,m))
      combine<CM>(rhs(i,j,k,l,m),lhs(i,j,k,l,m),alpha,beta);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    if (mask(i,j,k,l,m,n))
      combine<CM>(rhs(i,j,k,l,m,n),lhs(i,j,k,l,m,n),alpha,beta);
  }

  MaskView mask; 
  ST alpha;
  ST beta;
  LhsView lhs;
  RhsView rhs;
};

template<CombineMode CM, typename LhsView, typename RhsView, typename MaskView, typename ST>
void
cvmh (LhsView lhs, RhsView rhs,
      ST alpha, ST beta,
      const std::vector<int>& dims,
      MaskView mask)
{
  CombineViewsMaskedHelper <CM, LhsView, RhsView,  ST, MaskView> helper;
  helper.lhs = lhs;
  helper.rhs = rhs;
  helper.alpha = alpha;
  helper.beta = beta;
  helper.mask = mask;
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
      EKAT_ERROR_MSG ("Unsupported rank! Should be in [0,6].\n");
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if constexpr (N==0) {
      if (mask())
        lhs() = value;
    } else {
      if (mask(i))
        lhs(i) = value;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    if (mask(i,j))
      lhs(i,j) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    if (mask(i,j,k))
      lhs(i,j,k) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    if (mask(i,j,k,l))
      lhs(i,j,k,l) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m) const {
    if (mask(i,j,k,l,m))
      lhs(i,j,k,l,m) = value;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int m, int n) const {
    if (mask(i,j,k,l,m,n))
      lhs(i,j,k,l,m,n) = value;
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
     MaskView mask)
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
update_masked (const Field& x, const ST alpha, const ST beta, const Field& mask)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where one of the two views is strided (or both)
  const auto x_lr_ok = x.get_header().get_alloc_properties().allows_layout_right();
  const auto y_lr_ok = get_header().get_alloc_properties().allows_layout_right();
  const auto m_lr_ok = mask.get_header().get_alloc_properties().allows_layout_right();
  switch (layout.rank()) {
    case 0:
      if (x_lr_ok and y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST>(),x.get_view<const XST>(),
                                    alpha,beta,dims,mask.get_view<const int>())
                : details::cvmh<CM>(get_view<ST>(),x.get_view<const XST>(),
                                    alpha,beta,dims,mask.get_strided_view<const int>());
      else if (x_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST>(),x.get_view<const XST>(),
                                    alpha,beta,dims,mask.get_view<const int>())
                : details::cvmh<CM>(get_strided_view<ST>(),x.get_view<const XST>(),
                                    alpha,beta,dims,mask.get_strided_view<const int>());
      else if (y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST>(),x.get_strided_view<const XST>(),
                                    alpha,beta,dims,mask.get_view<const int>())
                : details::cvmh<CM>(get_view<ST>(),x.get_strided_view<const XST>(),
                                    alpha,beta,dims,mask.get_strided_view<const int>());
      else
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST>(),x.get_strided_view<const XST>(),
                                    alpha,beta,dims,mask.get_view<const int>())
                : details::cvmh<CM>(get_strided_view<ST>(),x.get_strided_view<const XST>(),
                                    alpha,beta,dims,mask.get_strided_view<const int>());
      break;
    case 1:
      if (x_lr_ok and y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST*>(),x.get_view<const XST*>(),
                                    alpha,beta,dims,mask.get_view<const int*>())
                : details::cvmh<CM>(get_view<ST*>(),x.get_view<const XST*>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*>());
      else if (x_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST*>(),x.get_view<const XST*>(),
                                    alpha,beta,dims,mask.get_view<const int*>())
                : details::cvmh<CM>(get_strided_view<ST*>(),x.get_view<const XST*>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*>());
      else if (y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST*>(),x.get_strided_view<const XST*>(),
                                    alpha,beta,dims,mask.get_view<const int*>())
                : details::cvmh<CM>(get_view<ST*>(),x.get_strided_view<const XST*>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*>());
      else
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST*>(),x.get_strided_view<const XST*>(),
                                    alpha,beta,dims,mask.get_view<const int*>())
                : details::cvmh<CM>(get_strided_view<ST*>(),x.get_strided_view<const XST*>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*>());
      break;
    case 2:
      if (x_lr_ok and y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST**>(),x.get_view<const XST**>(),
                                    alpha,beta,dims,mask.get_view<const int**>())
                : details::cvmh<CM>(get_view<ST**>(),x.get_view<const XST**>(),
                                    alpha,beta,dims,mask.get_strided_view<const int**>());
      else if (x_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST**>(),x.get_view<const XST**>(),
                                    alpha,beta,dims,mask.get_view<const int**>())
                : details::cvmh<CM>(get_strided_view<ST**>(),x.get_view<const XST**>(),
                                    alpha,beta,dims,mask.get_strided_view<const int**>());
      else if (y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST**>(),x.get_strided_view<const XST**>(),
                                    alpha,beta,dims,mask.get_view<const int**>())
                : details::cvmh<CM>(get_view<ST**>(),x.get_strided_view<const XST**>(),
                                    alpha,beta,dims,mask.get_strided_view<const int**>());
      else
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST**>(),x.get_strided_view<const XST**>(),
                                    alpha,beta,dims,mask.get_view<const int**>())
                : details::cvmh<CM>(get_strided_view<ST**>(),x.get_strided_view<const XST**>(),
                                    alpha,beta,dims,mask.get_strided_view<const int**>());
      break;
    case 3:
      if (x_lr_ok and y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST***>(),x.get_view<const XST***>(),
                                    alpha,beta,dims,mask.get_view<const int***>())
                : details::cvmh<CM>(get_view<ST***>(),x.get_view<const XST***>(),
                                    alpha,beta,dims,mask.get_strided_view<const int***>());
      else if (x_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST***>(),x.get_view<const XST***>(),
                                    alpha,beta,dims,mask.get_view<const int***>())
                : details::cvmh<CM>(get_strided_view<ST***>(),x.get_view<const XST***>(),
                                    alpha,beta,dims,mask.get_strided_view<const int***>());
      else if (y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST***>(),x.get_strided_view<const XST***>(),
                                    alpha,beta,dims,mask.get_view<const int***>())
                : details::cvmh<CM>(get_view<ST***>(),x.get_strided_view<const XST***>(),
                                    alpha,beta,dims,mask.get_strided_view<const int***>());
      else
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST***>(),x.get_strided_view<const XST***>(),
                                    alpha,beta,dims,mask.get_view<const int***>())
                : details::cvmh<CM>(get_strided_view<ST***>(),x.get_strided_view<const XST***>(),
                                    alpha,beta,dims,mask.get_strided_view<const int***>());
      break;
    case 4:
      if (x_lr_ok and y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST****>(),x.get_view<const XST****>(),
                                    alpha,beta,dims,mask.get_view<const int****>())
                : details::cvmh<CM>(get_view<ST****>(),x.get_view<const XST****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int****>());
      else if (x_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST****>(),x.get_view<const XST****>(),
                                    alpha,beta,dims,mask.get_view<const int****>())
                : details::cvmh<CM>(get_strided_view<ST****>(),x.get_view<const XST****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int****>());
      else if (y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST****>(),x.get_strided_view<const XST****>(),
                                    alpha,beta,dims,mask.get_view<const int****>())
                : details::cvmh<CM>(get_view<ST****>(),x.get_strided_view<const XST****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int****>());
      else
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST****>(),x.get_strided_view<const XST****>(),
                                    alpha,beta,dims,mask.get_view<const int****>())
                : details::cvmh<CM>(get_strided_view<ST****>(),x.get_strided_view<const XST****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int****>());
      break;
    case 5:
      if (x_lr_ok and y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST*****>(),x.get_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_view<const int*****>())
                : details::cvmh<CM>(get_view<ST*****>(),x.get_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*****>());
      else if (x_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST*****>(),x.get_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_view<const int*****>())
                : details::cvmh<CM>(get_strided_view<ST*****>(),x.get_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*****>());
      else if (y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST*****>(),x.get_strided_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_view<const int*****>())
                : details::cvmh<CM>(get_view<ST*****>(),x.get_strided_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*****>());
      else
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST*****>(),x.get_strided_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_view<const int*****>())
                : details::cvmh<CM>(get_strided_view<ST*****>(),x.get_strided_view<const XST*****>(),
                                    alpha,beta,dims,mask.get_strided_view<const int*****>());
      break;
    case 6:
      if (x_lr_ok and y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST******>(),x.get_view<const XST******>(),
                                    alpha,beta,dims,mask.get_view<const int******>())
                : details::cvmh<CM>(get_view<ST******>(),x.get_view<const XST******>(),
                                    alpha,beta,dims,mask.get_strided_view<const int******>());
      else if (x_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST******>(),x.get_view<const XST******>(),
                                    alpha,beta,dims,mask.get_view<const int******>())
                : details::cvmh<CM>(get_strided_view<ST******>(),x.get_view<const XST******>(),
                                    alpha,beta,dims,mask.get_strided_view<const int******>());
      else if (y_lr_ok)
        m_lr_ok ? details::cvmh<CM>(get_view<ST******>(),x.get_strided_view<const XST******>(),
                                    alpha,beta,dims,mask.get_view<const int******>())
                : details::cvmh<CM>(get_view<ST******>(),x.get_strided_view<const XST******>(),
                                    alpha,beta,dims,mask.get_strided_view<const int******>());
      else
        m_lr_ok ? details::cvmh<CM>(get_strided_view<ST******>(),x.get_strided_view<const XST******>(),
                                    alpha,beta,dims,mask.get_view<const int******>())
                : details::cvmh<CM>(get_strided_view<ST******>(),x.get_strided_view<const XST******>(),
                                    alpha,beta,dims,mask.get_strided_view<const int******>());
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_masked.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n"
          " - m name: " + mask.name() + "\n");
  }
  Kokkos::fence();
}

template<typename ST>
void Field::deep_copy_masked (const ST value, const Field& mask)
{
  const auto& layout = get_header().get_identifier().get_layout();
  const auto  rank   = layout.rank();
  const auto& dims   = layout.dims();
  const auto lr_ok = get_header().get_alloc_properties().allows_layout_right();
  const auto m_lr_ok = mask.get_header().get_alloc_properties().allows_layout_right();

  switch (rank) {
    case 0:
      if (lr_ok)
        return m_lr_ok ? details::svm(get_view<ST>(),value,dims,mask.get_view<const int>())
                       : details::svm(get_view<ST>(),value,dims,mask.get_strided_view<const int>());
      else
        return m_lr_ok ? details::svm(get_strided_view<ST>(),value,dims,mask.get_view<const int>())
                       : details::svm(get_strided_view<ST>(),value,dims,mask.get_strided_view<const int>());
      break;
    case 1:
      if (lr_ok)
        return m_lr_ok ? details::svm(get_view<ST*>(),value,dims,mask.get_view<const int*>())
                       : details::svm(get_view<ST*>(),value,dims,mask.get_strided_view<const int*>());
      else
        return m_lr_ok ? details::svm(get_strided_view<ST*>(),value,dims,mask.get_view<const int*>())
                       : details::svm(get_strided_view<ST*>(),value,dims,mask.get_strided_view<const int*>());
      break;
    case 2:
      if (lr_ok)
        return m_lr_ok ? details::svm(get_view<ST**>(),value,dims,mask.get_view<const int**>())
                       : details::svm(get_view<ST**>(),value,dims,mask.get_strided_view<const int**>());
      else
        return m_lr_ok ? details::svm(get_strided_view<ST**>(),value,dims,mask.get_view<const int**>())
                       : details::svm(get_strided_view<ST**>(),value,dims,mask.get_strided_view<const int**>());
      break;
    case 3:
      if (lr_ok)
        return m_lr_ok ? details::svm(get_view<ST***>(),value,dims,mask.get_view<const int***>())
                       : details::svm(get_view<ST***>(),value,dims,mask.get_strided_view<const int***>());
      else
        return m_lr_ok ? details::svm(get_strided_view<ST***>(),value,dims,mask.get_view<const int***>())
                       : details::svm(get_strided_view<ST***>(),value,dims,mask.get_strided_view<const int***>());
      break;
    case 4:
      if (lr_ok)
        return m_lr_ok ? details::svm(get_view<ST****>(),value,dims,mask.get_view<const int****>())
                       : details::svm(get_view<ST****>(),value,dims,mask.get_strided_view<const int****>());
      else
        return m_lr_ok ? details::svm(get_strided_view<ST****>(),value,dims,mask.get_view<const int****>())
                       : details::svm(get_strided_view<ST****>(),value,dims,mask.get_strided_view<const int****>());
      break;
    case 5:
      if (lr_ok)
        return m_lr_ok ? details::svm(get_view<ST*****>(),value,dims,mask.get_view<const int*****>())
                       : details::svm(get_view<ST*****>(),value,dims,mask.get_strided_view<const int*****>());
      else
        return m_lr_ok ? details::svm(get_strided_view<ST*****>(),value,dims,mask.get_view<const int*****>())
                       : details::svm(get_strided_view<ST*****>(),value,dims,mask.get_strided_view<const int*****>());
      break;
    case 6:
      if (lr_ok)
        return m_lr_ok ? details::svm(get_view<ST******>(),value,dims,mask.get_view<const int******>())
                       : details::svm(get_view<ST******>(),value,dims,mask.get_strided_view<const int******>());
      else
        return m_lr_ok ? details::svm(get_strided_view<ST******>(),value,dims,mask.get_view<const int******>())
                       : details::svm(get_strided_view<ST******>(),value,dims,mask.get_strided_view<const int******>());
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_MASKED_HPP
