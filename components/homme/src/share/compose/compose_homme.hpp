#ifndef INCLUDE_COMPOSE_HOMME_HPP
#define INCLUDE_COMPOSE_HOMME_HPP

#include "compose.hpp"
#include "compose_kokkos.hpp"

#include <cassert>
#include <memory>
#include <vector>

namespace homme {
typedef int Int;
typedef double Real;

namespace ko = Kokkos;

// Fortran array wrappers with Fortran index order.
template <typename T> using FA1 = ko::View<T*,     ko::LayoutLeft, ko::HostSpace>;
template <typename T> using FA2 = ko::View<T**,    ko::LayoutLeft, ko::HostSpace>;
template <typename T> using FA3 = ko::View<T***,   ko::LayoutLeft, ko::HostSpace>;
template <typename T> using FA4 = ko::View<T****,  ko::LayoutLeft, ko::HostSpace>;
template <typename T> using FA5 = ko::View<T*****, ko::LayoutLeft, ko::HostSpace>;

template <typename MT> using DepPoints =
  ko::View<Real****, ko::LayoutRight, typename MT::DDT>;
template <typename MT> using QExtrema =
  ko::View<Real****, ko::LayoutRight, typename MT::DDT>;
  
template <typename MT> using DepPointsH = typename DepPoints<MT>::HostMirror;
template <typename MT> using QExtremaH = typename QExtrema<MT>::HostMirror;

template <typename MT> using QExtremaHConst = ko::Const<QExtremaH<MT> >;
template <typename MT> using QExtremaConst = ko::Const<QExtrema<MT> >;

struct Cartesian3D { Real x, y, z; };

template <typename VT> COMPOSE_FORCEINLINE_FUNCTION
typename VT::value_type& idx_qext (const VT& qe, int ie, int q, int g, int lev) {
#ifdef COMPOSE_PORT
  return qe(ie,q,g,lev);
#else
  return qe(ie,q,lev,g);
#endif
}

template <typename T, int rank_>
struct HommeFormatSubArray {
  enum : bool { type_tag_HommeFormatSubArray = true };
  enum : int { rank = rank_ };
  typedef T value_type;

  HommeFormatSubArray (T* data_, Int np2_, Int nlev_ = -1, Int qsized_ = -1,
                       Int ntimelev_ = -1)
    : data(data_), nlev(nlev_), qsized(qsized_), ntimelev(ntimelev_)
  {
    assert(np2_ == np2);
  }

  COMPOSE_FORCEINLINE_FUNCTION
  T& operator() (const Int& i) const {
    static_assert(rank == 1, "rank 1 array");
    assert(i >= 0);
    return data[i];
  }
  COMPOSE_FORCEINLINE_FUNCTION 
  T& operator() (const Int& k, const Int& lev) const {
    static_assert(rank == 2, "rank 2 array");
    assert(k >= 0);
    assert(lev >= 0);
    check(k, lev);
    return data[lev*np2 + k];
  }
  COMPOSE_FORCEINLINE_FUNCTION 
  T& operator() (const Int& q_or_timelev, const Int& k, const Int& lev) const {
    static_assert(rank == 3, "rank 3 array");
    assert(q_or_timelev >= 0);
    assert(k >= 0);
    assert(lev >= 0);
    check(k, lev, q_or_timelev);
    return data[(q_or_timelev*nlev + lev)*np2 + k];
  }
  COMPOSE_FORCEINLINE_FUNCTION 
  T& operator() (const Int& timelev, const Int& q, const Int& k, const Int& lev) const {
    static_assert(rank == 4, "rank 4 array");
    assert(timelev >= 0);
    assert(q >= 0);
    assert(k >= 0);
    assert(lev >= 0);
    check(k, lev, q, timelev);
    return data[((timelev*qsized + q)*nlev + lev)*np2 + k];
  }

private:
  static const int np2 = 16;  
  T* data;
  const Int nlev, qsized, ntimelev;

  COMPOSE_FORCEINLINE_FUNCTION
  void check (Int k = -1, Int lev = -1, Int q_or_timelev = -1,
              Int timelev = -1) const {
#ifdef COMPOSE_BOUNDS_CHECK
    if (k >= 0) assert(k < np2);
    if (lev >= 0) assert(lev < nlev);
    if (q_or_timelev >= 0) {
      if (qsized < 0)
        assert(q_or_timelev < ntimelev);
      else
        assert(q_or_timelev < qsized);
    }
    if (timelev >= 0) assert(timelev < ntimelev);
#endif    
  }
};

template <typename T, int rank_>
struct HommeFormatArray {
  enum : bool { type_tag_HommeFormatArray = true };
  enum : int { rank = rank_ };
  typedef T value_type;
  typedef HommeFormatSubArray<T,rank_-1> subview_type;

  HommeFormatArray (Int nelemd, Int np2_, Int nlev_ = -1, Int qsized_ = -1,
                    Int ntimelev_ = -1)
    : nlev(nlev_), qsized(qsized_), ntimelev(ntimelev_)
  {
    assert(np2_ == np2);
    ie_data_ptr.resize(nelemd);
  }

  void set_ie_ptr (const Int ie, T* ptr) {
    check(ie);
    ie_data_ptr[ie] = ptr;
  }

  COMPOSE_FORCEINLINE_FUNCTION
  subview_type get_sub (const Int& ie) const {
    return subview_type(ie_data_ptr[ie], np2, nlev, qsized, ntimelev);
  }

  COMPOSE_FORCEINLINE_FUNCTION
  T& operator() (const Int& ie, const Int& i) const {
    static_assert(rank == 2, "rank 2 array");
    // These routines are not used on the GPU, but they can be called from
    // KOKKOS_FUNCTIONs on CPU in GPU builds. Avoid nvcc warnings as follows:
#if defined __CUDA_ARCH__ || defined __HIP_DEVICE_COMPILE__
    return unused();
#else
    assert(i >= 0);
    assert(ie_data_ptr[ie]);
    return *(ie_data_ptr[ie] + i);
#endif
  }
  COMPOSE_FORCEINLINE_FUNCTION 
  T& operator() (const Int& ie, const Int& k, const Int& lev) const {
    static_assert(rank == 3, "rank 3 array");
#if defined __CUDA_ARCH__ || defined __HIP_DEVICE_COMPILE__
    return unused();
#else
    assert(k >= 0);
    assert(lev >= 0);
    assert(ie_data_ptr[ie]);
    check(ie, k, lev);
    return *(ie_data_ptr[ie] + lev*np2 + k);
#endif
  }
  COMPOSE_FORCEINLINE_FUNCTION 
  T& operator() (const Int& ie, const Int& q_or_timelev, const Int& k,
                 const Int& lev) const {
    static_assert(rank == 4, "rank 4 array");
#if defined __CUDA_ARCH__ || defined __HIP_DEVICE_COMPILE__
    return unused();
#else
    assert(q_or_timelev >= 0);
    assert(k >= 0);
    assert(lev >= 0);
    assert(ie_data_ptr[ie]);
    check(ie, k, lev, q_or_timelev);
    return *(ie_data_ptr[ie] + (q_or_timelev*nlev + lev)*np2 + k);
#endif
  }
  COMPOSE_FORCEINLINE_FUNCTION 
  T& operator() (const Int& ie, const Int& timelev, const Int& q, const Int& k,
                 const Int& lev) const {
    static_assert(rank == 5, "rank 4 array");
#if defined __CUDA_ARCH__ || defined __HIP_DEVICE_COMPILE__
    return unused();
#else
    assert(timelev >= 0);
    assert(q >= 0);
    assert(k >= 0);
    assert(lev >= 0);
    assert(ie_data_ptr[ie]);
    check(ie, k, lev, q, timelev);
    return *(ie_data_ptr[ie] + ((timelev*qsized + q)*nlev + lev)*np2 + k);
#endif
  }

private:
  static const int np2 = 16;
  std::vector<T*> ie_data_ptr;
  const Int nlev, qsized, ntimelev;

#ifdef COMPOSE_ENABLE_GPU
  COMPOSE_INLINE_FUNCTION static T& unused () {
    static T unused = 0;
    assert(0);
    return unused;
  }
#endif

  COMPOSE_FORCEINLINE_FUNCTION
  void check (Int ie, Int k = -1, Int lev = -1, Int q_or_timelev = -1,
              Int timelev = -1) const {
#if defined COMPOSE_BOUNDS_CHECK && ! (defined __CUDA_ARCH__ || defined __HIP_DEVICE_COMPILE__)
    assert(ie >= 0 && ie < static_cast<Int>(ie_data_ptr.size()));
    if (k >= 0) assert(k < np2);
    if (lev >= 0) assert(lev < nlev);
    if (q_or_timelev >= 0) {
      if (qsized < 0)
        assert(q_or_timelev < ntimelev);
      else
        assert(q_or_timelev < qsized);
    }
    if (timelev >= 0) assert(timelev < ntimelev);
#endif    
  }
};

template <typename HFA> inline HFA& unmanaged (
  HFA& s, typename std::enable_if<HFA::type_tag_HommeFormatArray>::type* = 0)
{ return s; }

template <typename HFA> inline typename HFA::subview_type subview_ie (
  const Int ie, HFA& s, typename std::enable_if<HFA::type_tag_HommeFormatArray>::type* = 0)
{ return s.get_sub(ie); }

// Qdp, dp, Q
template <typename MT>
struct TracerArrays {
  typedef std::shared_ptr<TracerArrays<MT> > Ptr;

  const Int nelemd, nlev, np, np2, qsize, qsized;
  HommeFormatArray<const Real,2> pspheremp;
  HommeFormatArray<const Real,3> pdp;
  HommeFormatArray<const Real,4> pdp3d;
  HommeFormatArray<Real,5> pqdp;
  HommeFormatArray<Real,4> pq;
  Int n0_qdp, n1_qdp, np1;

#if defined COMPOSE_PORT
  template <typename Datatype> using View =
    ko::View<Datatype, ko::LayoutRight, typename MT::DDT>;
  View<Real**> spheremp;
  View<Real***>   dp;   // elem%derived%dp
  View<Real****>  dp3d; // elem%state%dp3d or the sl3d equivalent
  View<Real*****> qdp;  // elem%state%Qdp(:,:,:,:,:)
  View<Real****>  q;    // elem%state%Q
  DepPoints<MT> dep_points, vnode, vdep;
  QExtrema<MT> q_min, q_max;
  void alloc_if_not();
#else
  HommeFormatArray<const Real,2>& spheremp;
  HommeFormatArray<const Real,3>& dp;
  HommeFormatArray<const Real,4>& dp3d;
  HommeFormatArray<Real,5>& qdp;
  HommeFormatArray<Real,4>& q;
#endif

  TracerArrays(Int nelemd, Int nlev, Int np, Int qsize, Int qsized);
  TracerArrays(const TracerArrays<MT>&) = delete;
  TracerArrays& operator=(const TracerArrays<MT>&) = delete;
};

template <typename DataType>
using TracerView = ko::View<DataType, ko::LayoutRight, ko::MachineTraits::DDT>;
template <typename T> COMPOSE_INLINE_FUNCTION TracerView<T*>
subview_ie (const Int ie, const TracerView<T**>& s)
{ return TracerView<T*>(&s(ie,0), s.extent(1)); }
template <typename T> COMPOSE_INLINE_FUNCTION TracerView<T**>
subview_ie (const Int ie, const TracerView<T***>& s)
{ return TracerView<T**>(&s(ie,0,0), s.extent(1), s.extent(2)); }
template <typename T> COMPOSE_INLINE_FUNCTION TracerView<T***>
subview_ie (const Int ie, const TracerView<T****>& s)
{ return TracerView<T***>(&s(ie,0,0,0), s.extent(1), s.extent(2), s.extent(3)); }
template <typename T> COMPOSE_INLINE_FUNCTION TracerView<T****>
subview_ie (const Int ie, const TracerView<T*****>& s)
{ return TracerView<T****>(&s(ie,0,0,0,0), s.extent(1), s.extent(2), s.extent(3), s.extent(4)); }

template <typename MT>
void sl_traj_h2d(TracerArrays<MT>& ta, Real* dep_points, Real* vnode, Real* vdep,
                 Int ndim);

template <typename MT>
void sl_traj_d2h(const TracerArrays<MT>& ta, Real* dep_points, Real* vnode,
                 Real* vdep, Int ndim);

template <typename MT>
void sl_h2d(TracerArrays<MT>& ta, bool transfer, Real* dep_points, Int ndim);

template <typename MT>
void sl_d2h(const TracerArrays<MT>& ta, bool transfer, Real* dep_points, Int ndim,
            Real* minq, Real* maxq);

template <typename MT>
void cedr_h2d(const TracerArrays<MT>& ta, bool transfer);

template <typename MT>
void cedr_d2h(const TracerArrays<MT>& ta, bool transfer);

TracerArrays<ko::MachineTraits>::Ptr
init_tracer_arrays(Int nelemd, Int nlev, Int np, Int qsize, Int qsize_d);

TracerArrays<ko::MachineTraits>::Ptr get_tracer_arrays();

void delete_tracer_arrays();

} // namespace homme

#endif
