/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#ifndef HOMMEXX_COMPOSE_TRANSPORT_IMPL_ENHANCED_TRAJECTORY_IMPL_HPP
#define HOMMEXX_COMPOSE_TRANSPORT_IMPL_ENHANCED_TRAJECTORY_IMPL_HPP

#include "ComposeTransportImpl.hpp"

#include "compose_hommexx.hpp"

namespace Homme {

// Organize and interpolate velocity snapshots from the dynamics. The result is
// a set of snapshots to be used in computing the SL trajectories. To
// distinguish between these types of snapshots, we refer to "dynamics
// snapshots" (available at vertical remap time steps) and "internal snapshots"
// (interpolated to tracer trajectory substeps).
//
// Parameter short names:
//   dtf = dt_tracer_factor
//   drf = dt_remap_factor
//   nsub = semi_lagrange_trajectory_nsubstep
//   nvel = semi_lagrange_trajectory_nvelocity
struct ComposeTransportImpl::VelocityRecord {
  VelocityRecord () {}
  VelocityRecord (const int dtf, const int drf, const int nsub, const int nvel)
  { init(dtf, drf, nsub, nvel); }

  void init(const int dtf, const int drf, const int nsub, const int nvel);

  int dtf () const { return _dtf; }
  int drf () const { return _drf; }
  int nsub () const { return _nsub; }
  int nvel () const { return _nvel; }

  // Times to which velocity slots i in 0:nvel-1 correspond, in reference time
  // [0,dtf].
  Real t_vel(const int i) const;

  // For n = 0:dtf-1, obs_slots(n,0:1) = [slot1, slot2], -1 if unused. These are
  // the slots to which velocity sample n contributes. obs_slots(dtf-1,:) is
  // always -1.
  int obs_slots(const int n, const int k) const;

  // For n = 0:dtf-1, obs_wts(n,0:1) = [wt1, wt2], 0 if unused.
  Real obs_wts(const int n, const int k) const;

  // Substep end point i in 0:nsub uses velocity slots run_step(i),
  // run_step(i)-1.
  int run_step(const int i) const;

private:
  int _dtf = -1, _drf = -1, _nsub = -1, _nvel = -1;
  std::vector<int> _obs_slots, _run_step;
  std::vector<Real> _t_vel, _obs_wts;

  Real& t_vel(const int i);
  int& obs_slots(const int n, const int k);
  Real& obs_wts(const int n, const int k);
  int& run_step(const int i);
};

namespace { // anon

using cti = ComposeTransportImpl;
using CTI = ComposeTransportImpl;
using CSelNlev  = cti::CSNlev;
using CRelNlev  = cti::CRNlev;
using CSelNlevp = cti::CSNlevp;
using CRelNlevp = cti::CRNlevp;
using CS2elNlev = cti::CS2Nlev;
using CR2elNlev = cti::CR2Nlev;
using SelNlev   = cti::SNlev;
using RelNlev   = cti::RNlev;
using SelNlevp  = cti::SNlevp;
using RelNlevp  = cti::RNlevp;
using S2elNlev  = cti::S2Nlev;
using R2elNlev  = cti::R2Nlev;
using S2elNlevp = cti::S2Nlevp;

using RelV = ExecViewUnmanaged<Real[NP][NP]>;
using CRelV = typename ViewConst<RelV>::type;

template <int N> using SelNV = ExecViewUnmanaged<Scalar[NP][NP][N]>;
template <int N> using CSelNV = typename ViewConst<SelNV<N>>::type;

template <int N> using RelNV = ExecViewUnmanaged<Real[NP][NP][N]>;
template <int N> using CRelNV = typename ViewConst<RelNV<N>>::type;

template <int N> using RNV = ExecViewUnmanaged<Real[N]>;
template <int N> using CRNV = typename ViewConst<RNV<N>>::type;
using RNlevp = RNV<cti::num_phys_lev+1>;
using CRNlevp = CRNV<cti::num_phys_lev+1>;

using RnV = ExecViewUnmanaged<Real*>;
using CRnV = ExecViewUnmanaged<const Real*>;
using SnV = ExecViewUnmanaged<Scalar*>;
using CSnV = ExecViewUnmanaged<const Scalar*>;

template <int N> using SNV = ExecViewUnmanaged<Scalar[N]>;
template <int N> using CSNV = typename ViewConst<SNV<N>>::type;

using RelnV = ExecViewUnmanaged<Real***>;
using CRelnV = ExecViewUnmanaged<const Real***>;
using SelnV = ExecViewUnmanaged<Scalar***>;
using CSelnV = ExecViewUnmanaged<const Scalar***>;

// Helper functions to move between various array data structures and assert
// things about them.

KOKKOS_INLINE_FUNCTION
static int calc_npack (const int nscal) {
  return (nscal + cti::packn - 1) / VECTOR_SIZE;
}

KOKKOS_INLINE_FUNCTION
static int calc_nscal (const int npack) {
  return npack * VECTOR_SIZE;
}

KOKKOS_INLINE_FUNCTION
RnV getcol (const RelnV& a, const int i, const int j) {
  return Kokkos::subview(a,i,j,Kokkos::ALL);
}

KOKKOS_INLINE_FUNCTION
CRnV getcolc (const CRelnV& a, const int i, const int j) {
  return Kokkos::subview(a,i,j,Kokkos::ALL);
}

KOKKOS_INLINE_FUNCTION
RelnV elp2r (const SelnV& p) {
  return RelnV(cti::pack2real(p), NP, NP, calc_nscal(p.extent_int(2)));
}

KOKKOS_INLINE_FUNCTION
CRelnV elp2r (const CSelnV& p) {
  return CRelnV(cti::cpack2real(p), NP, NP, calc_nscal(p.extent_int(2)));
}

KOKKOS_INLINE_FUNCTION
RelnV p2rel (Scalar* data, const int nlev) {
  return RelnV(cti::pack2real(data), NP, NP, nlev);
}

KOKKOS_INLINE_FUNCTION
void assert_eln (const CRelnV& a, const int nlev) {
  assert(a.extent_int(0) >= NP);
  assert(a.extent_int(1) >= NP);
  assert(a.extent_int(2) >= nlev);
}

KOKKOS_INLINE_FUNCTION
void assert_eln (const CSelnV& a, const int nlev) {
  assert(a.extent_int(0) >= NP);
  assert(a.extent_int(1) >= NP);
  assert(calc_nscal(a.extent_int(2)) >= nlev);
}

// Structs to manage access to internal velocity snapshots at the end points of
// a substep interval.
//   dp1,2 and v1,2 are on Eulerian levels. dp,v1 is from time t1 < t2.

using  DpSnap   = ExecViewUnmanaged<      Scalar**    [NP][NP][NUM_LEV]>;
using   VSnap   = ExecViewUnmanaged<      Scalar** [2][NP][NP][NUM_LEV]>;
using CDpSnap   = ExecViewUnmanaged<const Scalar**    [NP][NP][NUM_LEV]>;
using  CVSnap   = ExecViewUnmanaged<const Scalar** [2][NP][NP][NUM_LEV]>;
using CDpSnapEl = ExecViewUnmanaged<const Scalar      [NP][NP][NUM_LEV]>;
using  CVSnapEl = ExecViewUnmanaged<const Scalar   [2][NP][NP][NUM_LEV]>;

// This is the simple case, for nvelocity = 2. We have velocity snapshots only
// at the tracer time step end points; dynamics and internal are the same. We
// linearly interpolate them to give values anywhere within the time step. In
// particular, for a particular substep's interval, we linearly interpolate
// using alpha[t] for t = 0,1 the start and end of the substep interval.
struct EndpointSnapshots {
  const Real alpha[2];
  const int idxs[2];
  const CDpSnap dps[2];
  const CVSnap vs[2];

  // Use subview access for efficiency.
  struct Element {
    const Real alpha[2];
    const CDpSnapEl dps[2];
    const CVSnapEl vs[2];

    KOKKOS_INLINE_FUNCTION
    Element (const EndpointSnapshots& s, const int ie)
      : alpha{s.alpha[0], s.alpha[1]},
        dps{Homme::subview(s.dps[0], ie, s.idxs[0]), Homme::subview(s.dps[1], ie, s.idxs[1])},
        vs{Homme::subview(s.vs[0], ie, s.idxs[0]), Homme::subview(s.vs[1], ie, s.idxs[1])}
    {}

    // Direct access.

    KOKKOS_INLINE_FUNCTION
    Real get_dp_real (const int t, const int i, const int j, const int k) const {
      return dps[t](i,j, k / VECTOR_SIZE)[k % VECTOR_SIZE];
    }

    // Linear combinations.

    KOKKOS_INLINE_FUNCTION
    Scalar combine_v (const int t, const int d, const int i, const int j, const int k) const {
      return (1 - alpha[t])*vs[0](d,i,j,k) + alpha[t]*vs[1](d,i,j,k);
    }

    KOKKOS_INLINE_FUNCTION
    Scalar combine_vdp (const int t, const int d, const int i, const int j, const int k) const {
      return ((1 - alpha[t])*vs[0](d,i,j,k)*dps[0](i,j,k)
              +    alpha[t] *vs[1](d,i,j,k)*dps[1](i,j,k));
    }
  };

  EndpointSnapshots (const Real alpha_[2],
                     const CDpSnap& dp1, const CVSnap& v1, const int idx1,
                     const CDpSnap& dp2, const CVSnap& v2, const int idx2)
    : alpha{alpha_[0], alpha_[1]}, idxs{idx1, idx2}, dps{dp1, dp2}, vs{v1, v2}
  {}

  KOKKOS_INLINE_FUNCTION Real get_alpha (const int t) const { return alpha[t]; }

  KOKKOS_INLINE_FUNCTION Element get_element (const int ie) const {
    return Element(*this, ie);
  }
};

// This is the more complicated case: we have dynamics velocity snapshots at the
// time step end points, plus internal snapshots that may be interpolated from
// multiple dynamics snapshots. alpha is no longer needed; values are trivially
// set to 0 or 1. Interpolation of snapshots here takes the place of alpha.
struct ManySnapshots {
  Real beta[2];
  int idxs[4];
  CDpSnap dps[4];
  CVSnap vs[4];

  struct Element {
    const Real beta[2];
    const CDpSnapEl dps[4];
    const CVSnapEl vs[4];

    KOKKOS_INLINE_FUNCTION
    Element (const ManySnapshots& s, const int ie)
      : beta{s.beta[0], s.beta[1]},
        dps{Homme::subview(s.dps[0], ie, s.idxs[0]), Homme::subview(s.dps[1], ie, s.idxs[1]),
            Homme::subview(s.dps[2], ie, s.idxs[2]), Homme::subview(s.dps[3], ie, s.idxs[3])},
        vs{Homme::subview(s.vs[0], ie, s.idxs[0]), Homme::subview(s.vs[1], ie, s.idxs[1]),
           Homme::subview(s.vs[2], ie, s.idxs[2]), Homme::subview(s.vs[3], ie, s.idxs[3])}
    {}

    KOKKOS_INLINE_FUNCTION
    Real get_dp_real (const int t, const int i, const int j, const int k) const {
      const int os = 2*t, kp = k / VECTOR_SIZE, ks = k % VECTOR_SIZE;
      return ((1 - beta[t])*dps[os  ](i,j,kp)[ks]
              +    beta[t] *dps[os+1](i,j,kp)[ks]);
    }

    KOKKOS_INLINE_FUNCTION
    Scalar combine_v (const int t, const int d, const int i, const int j, const int k) const {
      const int os = 2*t;
      return ((1 - beta[t])*vs[os  ](d,i,j,k)
              +    beta[t] *vs[os+1](d,i,j,k));
    }

    KOKKOS_INLINE_FUNCTION
    Scalar combine_vdp (const int t, const int d, const int i, const int j, const int k) const {
      const int os = 2*t;
      return ((1 - beta[t])*vs[os  ](d,i,j,k)*dps[os  ](i,j,k)
              +    beta[t] *vs[os+1](d,i,j,k)*dps[os+1](i,j,k));
    }
  };

  ManySnapshots (const CTI::VelocityRecord& vrec,
                 // endpoint data
                 const CDpSnap& dp1, const CVSnap& v1, const int idx1,
                 const CDpSnap& dp2, const CVSnap& v2, const int idx2,
                 // interior data
                 const CDpSnap& dps_int, const CVSnap& vs_int,
                 // current substep
                 const int nsubstep, const int step) {
    const int end = vrec.nvel() - 1;
    for (int t = 0; t < 2; ++t) {
      const int substep_idx = nsubstep - (step+1) + t;
      const Real time = (substep_idx * vrec.t_vel(end))/nsubstep;
      const int k = vrec.run_step(step);
      const int os = 2*t;
      // Get the two relevant dynamics snapshots.
      idxs[os]   = k == 1   ? idx1 : (k-2);
      dps [os]   = k == 1   ? dp1  : dps_int;
      vs  [os]   = k == 1   ? v1   : vs_int;
      idxs[os+1] = k == end ? idx2 : (k-1);
      dps[os+1]  = k == end ? dp2  : dps_int;
      vs[os+1]   = k == end ? v2   : vs_int;
      assert(idxs[os] >= 0 and idxs[os+1] >= 0);
      assert(idxs[os] < dps[os].extent_int(1) and idxs[os+1] < dps[os+1].extent_int(1));
      // Parameter for the linear combination of the two dynamics snapshots to
      // make an internal snapshot.
      beta[t] = (time - vrec.t_vel(k-1))/(vrec.t_vel(k) - vrec.t_vel(k-1));
    }
  }

  KOKKOS_INLINE_FUNCTION Real get_alpha (const int t) const {
    return t == 0 ? 0 : 1;
  }

  KOKKOS_INLINE_FUNCTION Element get_element (const int ie) const {
    return Element(*this, ie);
  }
};

// Routines to interpolate in the vertical direction.

// Find the support for the linear interpolant.
//   For sorted ascending x[0:n] and x in [x[0], x[n-1]] with hint xi_idx,
// return i such that x[i] <= xi <= x[i+1].
//   This function is meant for the case that x_idx is very close to the
// support. If that isn't true, then this method is inefficient; binary search
// should be used instead.
template <typename ConstRealArray>
KOKKOS_FUNCTION static
int find_support (const int n, const ConstRealArray& x, const int x_idx,
                  const Real xi) {
  assert(xi >= x[0] and xi <= x[n-1]);
  // Handle the most common case.
  if (x_idx < n-1 and xi >= x[x_idx  ] and xi <= x[x_idx+1]) return x_idx;
  if (x_idx > 0   and xi >= x[x_idx-1] and xi <= x[x_idx  ]) return x_idx-1;
  // Move on to less common ones.
  const int max_step = max(x_idx, n-1 - x_idx);
  for (int step = 1; step <= max_step; ++step) {
    if (x_idx < n-1-step and xi >= x[x_idx+step  ] and xi <= x[x_idx+step+1])
      return x_idx+step;
    if (x_idx > step     and xi >= x[x_idx-step-1] and xi <= x[x_idx-step  ])
      return x_idx-step-1;
  }
  assert(false);
  return -1;
}

// Linear interpolation core formula.
template <typename XT, typename YT>
KOKKOS_FUNCTION Real
linterp (const int n, const XT& x, const YT& y, const int x_idx, const Real xi) {
  const auto isup = find_support(n, x, x_idx, xi);
  const Real a = (xi - x[isup])/(x[isup+1] - x[isup]);
  return (1-a)*y[isup] + a*y[isup+1];
}

// Linear interpolation at the lowest level of team ||ism.
//   Range provides this ||ism over index 0 <= k < ni.
//   Interpolate y(x) to yi(xi).
//   x_idx_offset is added to k in the call to find_support.
//   Arrays should all have rank 1.
//   Notation: yi = I[y(x)](xi) is an interpolant constructed from y(x) and
//             evaluated at xi.
template <typename Range, typename XT, typename YT, typename XIT, typename YIT>
KOKKOS_FUNCTION void
linterp (const Range& range,
         const int n , const XT& x , const YT& y,
         const int ni, const XIT& xi, const YIT& yi,
         const int x_idx_offset = 0, const char* const caller = nullptr) {
#ifndef NDEBUG
  if (xi[0] < x[0] or xi[ni-1] > x[n-1]) {
    Kokkos::printf("linterp: xi out of bounds: %s %1.15e %1.15e %1.15e %1.15e\n",
                   caller ? caller : "NONE", x[0], xi[0], xi[ni-1], x[n-1]);
    assert(false);
  }
#endif
  assert(range.start ==  0);
  assert(range.end   == ni);
  const auto f = [&] (const int k) {
    yi[k] = linterp(n, x, y, k + x_idx_offset, xi[k]);
  };
  Kokkos::parallel_for(range, f);
}

// Compute Lagrangian levels at t1 on arrival column:
//     yi(i_os:) = I[y([eta(0),x,eta(1)])](xi(i_os:)),
// where both x and y are eta values located at midpoints or interfaces and on
// the reference or departure grids.
KOKKOS_FUNCTION void
eta_interp_eta (const KernelVariables& kv, const int nlev, const CRnV& hy_etai,
                const int n, const int os, const CRelnV& x, const CRnV& y,
                const RelnV& xwrk, const RnV& ywrk,
                // Use xi(i_os:), yi(i,j,i_os:).
                const int ni, const int i_os, const CRnV& xi, const RelnV& yi) {
  const auto& xbdy = xwrk;
  const auto& ybdy = ywrk;
  assert(hy_etai.extent_int(0) >= nlev+1);
  assert_eln(x, os + n);
  assert(y.extent_int(0) >= os + n);
  assert_eln(xbdy, n+2);
  assert(ybdy.extent_int(0) >= n+2);
  assert(xi.extent_int(0) >= i_os + ni);
  assert_eln(yi, i_os + ni);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr_ni = Kokkos::ThreadVectorRange(kv.team, ni);
  const auto tvr_np2 = Kokkos::ThreadVectorRange(kv.team, n+2);
  const auto f_y = [&] (const int k) {
    ybdy(k) = (k == 0   ? hy_etai(0) :
               k == n+1 ? hy_etai(nlev) :
               /**/       y(os+k-1));
  };
  Kokkos::parallel_for(Kokkos::TeamVectorRange(kv.team, n+2), f_y);
  kv.team_barrier();
  const auto f_x = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto g = [&] (const int k) {
      xbdy(i,j,k) = (k == 0   ? hy_etai(0) :
                     k == n+1 ? hy_etai(nlev) :
                     /**/       x(i,j,os+k-1));
    };
    Kokkos::parallel_for(tvr_np2, g);
  };
  Kokkos::parallel_for(ttr, f_x);
  kv.team_barrier();
  const auto f_linterp = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    linterp(tvr_ni,
            n+2, getcolc(xbdy,i,j), ybdy,
            ni, xi.data() + i_os, getcol(yi,i,j).data() + i_os,
            1, "eta_interp_eta");
  };
  Kokkos::parallel_for(ttr, f_linterp);
}

// Compute departure horizontal points corresponding to arrival Lagrangian level
// midpoints:
//     p_dep_mid(eta_arr_mid) = I[p_dep_mid(eta_ref_mid)](eta_arr_mid)
KOKKOS_FUNCTION void
eta_interp_horiz (const KernelVariables& kv, const int nlev,
                  const CRnV& hy_etai, const CRnV& x, const CRelnV& y,
                  const RnV& xwrk, const RelnV& ywrk,
                  const CRelnV& xi, const RelnV& yi) {
  const auto& xbdy = xwrk;
  const auto& ybdy = ywrk;
  assert(hy_etai.extent_int(0) >= nlev+1);
  assert(x.extent_int(0) >= nlev);
  assert_eln(y, nlev);
  assert(xbdy.extent_int(0) >= nlev+2);
  assert_eln(ybdy, nlev+2);
  assert_eln(xi, nlev);
  assert_eln(yi, nlev);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr_nlev = Kokkos::ThreadVectorRange(kv.team, nlev);
  const auto tvr_nlevp2 = Kokkos::ThreadVectorRange(kv.team, nlev+2);
  const auto f_x = [&] (const int k) {
    xbdy(k) = (k == 0      ? hy_etai(0) :
               k == nlev+1 ? hy_etai(nlev) :
               /**/          x(k-1));
  };
  Kokkos::parallel_for(Kokkos::TeamVectorRange(kv.team, nlev+2), f_x);
  kv.team_barrier();
  const auto f_y = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto g = [&] (const int k) {
      // Constant interp outside of the etam support.
      ybdy(i,j,k) = (k == 0      ? y(i,j,0) :
                     k == nlev+1 ? y(i,j,nlev-1) :
                     /**/          y(i,j,k-1));
    };
    Kokkos::parallel_for(tvr_nlevp2, g);
  };
  Kokkos::parallel_for(ttr, f_y);
  kv.team_barrier();
  const auto f_linterp = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    linterp(tvr_nlev,
            nlev+2, xbdy,            getcolc(ybdy,i,j),
            nlev,   getcolc(xi,i,j), getcol(yi,i,j),
            1, "eta_interp_horiz");
  };
  Kokkos::parallel_for(ttr, f_linterp);
}

/* Compute level pressure thickness given eta at interfaces using the following
   approximation:
          e = A(e) + B(e)
       p(e) = A(e) p0 + B(e) ps
            = e p0 + B(e) (ps - p0)
           a= e p0 + I[Bi(eref)](e) (ps - p0).
   Then dp = diff(p).
*/
KOKKOS_FUNCTION void
eta_to_dp (const KernelVariables& kv, const int nlev,
           const Real hy_ps0, const CRnV& hy_bi, const CRnV& hy_etai,
           const CRelV& ps, const CRelnV& etai, const RelnV& wrk,
           const RelnV& dp) {
  const int nlevp = nlev + 1;
  assert(hy_bi.extent_int(0) >= nlevp);
  assert(hy_etai.extent_int(0) >= nlevp);
  assert_eln(etai, nlevp);
  assert_eln(wrk, nlevp);
  assert_eln(dp, nlev);
  const auto& bi = wrk;
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr_linterp = Kokkos::ThreadVectorRange(kv.team, nlevp);
  const auto f_linterp = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    linterp(tvr_linterp,
            nlevp, hy_etai, hy_bi,
            nlevp, getcolc(etai,i,j), getcol(bi,i,j),
            0, "eta_to_dp");
  };
  Kokkos::parallel_for(ttr, f_linterp);
  kv.team_barrier();
  const auto tvr = Kokkos::ThreadVectorRange(kv.team, nlev);
  const auto f = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto dps = ps(i,j) - hy_ps0;
    const auto g = [&] (const int k) {
      dp(i,j,k) = ((etai(i,j,k+1) - etai(i,j,k))*hy_ps0 +
                   (bi(i,j,k+1) - bi(i,j,k))*dps);
    };
    Kokkos::parallel_for(tvr, g);
  };
  Kokkos::parallel_for(ttr, f);
}

/* Limit eta levels so their thicknesses, deta, are bounded below by 'low'.

   This method pulls mass only from intervals k that are larger than their
   reference value (deta(k) > deta_ref(k)), and only down to their reference
   value. This concentrates changes to intervals that, by having a lot more mass
   than usual, drive other levels negative, leaving all the other intervals
   unchanged.

   This selective use of mass provides enough to fulfill the needed mass.
   Inputs:
       m (deta): input mass
       r (deta_ref): level mass reference.
   Preconditions:
       (1) 0 <= low <= min r(i)
       (2) 1 = sum r(i) = sum(m(i)).
   Rewrite (2) as
       1 = sum_{m(i) >= r(i)} m(i) + sum_{m(i) < r(i)} m(i)
   and, thus,
       0 = sum_{m(i) >= r(i)} (m(i) - r(i)) + sum_{m(i) < r(i)} (m(i) - r(i)).
   Then
       sum_{m(i) >= r(i)} (m(i) - r(i))         (available mass to redistribute)
           = -sum_{m(i) < r(i)} (m(i) - r(i))
          >= -sum_{m(i) < lo  } (m(i) - r(i))
          >= -sum_{m(i) < lo  } (m(i) - lo  )   (mass to fill in).
   Thus, if the preconditions hold, then there's enough mass to redistribute.
 */
template <typename Range>
KOKKOS_FUNCTION void
deta_caas (const KernelVariables& kv, const Range& tvr,
           const CRnV& deta_ref, const Real low, const RnV& w,
           const RnV& deta) {
  const auto g1 = [&] (const int k, Kokkos::Real2& sums) {
    Real wk;
    if (deta(k) < low) {
      sums.v[0] += deta(k) - low;
      deta(k) = low;
      wk = 0;
    } else {
      wk = (deta(k) > deta_ref(k) ?
            deta(k) - deta_ref(k) :
            0);
    }
    sums.v[1] += wk;
    w(k) = wk;
  };
  Kokkos::Real2 sums;
  Dispatch<>::parallel_reduce(kv.team, tvr, g1, sums);
  const Real wneeded = sums.v[0];
  if (wneeded == 0) return;
  // Remove what is needed from the donors.
  const Real wavail = sums.v[1];
  const auto g2 = [&] (const int k) {
    deta(k) += wneeded*(w(k)/wavail);
  };
  Kokkos::parallel_for(tvr, g2);
}

// Wrapper to deta_caas. On input and output, eta contains the interface eta
// values, excluding the last one. On output, deta_caas has been applied, if
// necessary, to diff(eta(i,j,:)).
KOKKOS_FUNCTION int
limit_etai (const KernelVariables& kv, const int nlev, const CRnV& hy_etai,
            const CRnV& deta_ref, const Real deta_tol, const RelnV& wrk1,
            const RelnV& wrk2, const RelnV& eta) {
  assert(hy_etai.extent_int(0) >= nlev+1);
  assert(deta_ref.extent_int(0) >= nlev);
  const auto deta = wrk2;
  assert_eln(wrk1, nlev);
  assert_eln(deta, nlev);
  assert_eln(eta , nlev);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr = Kokkos::ThreadVectorRange(kv.team, nlev);
  // eta -> deta; limit deta if needed.
  const auto f1 = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto  etaij = getcolc( eta,i,j);
    const auto detaij = getcol (deta,i,j);
    const auto g1 = [&] (const int k, int& nbad) {
      const auto d = (k+1 == nlev ? hy_etai(nlev) - etaij(nlev-1) :
                      /**/          etaij(k+1) - etaij(k));
      const bool ok = d >= deta_tol;
      if (not ok) ++nbad;
      detaij(k) = d;
    };
    int nbad = 0;
    Dispatch<>::parallel_reduce(kv.team, tvr, g1, nbad);
    if (nbad == 0) {
      // Signal this column is fine.
      Kokkos::single(Kokkos::PerThread(kv.team), [&] () { detaij(0) = -1; });
      return;
    }
    deta_caas(kv, tvr, deta_ref, deta_tol, getcol(wrk1,i,j), detaij);
  };
  Kokkos::parallel_for(ttr, f1);
  kv.team_barrier();
  // deta -> eta; ignore columns where limiting wasn't needed.
  const auto f2 = [&] (const int idx, int& cnt) {
    const int i = idx / NP, j = idx % NP;
    const auto detaij = getcolc(deta,i,j);
    if (detaij(0) == -1) return;
    ++cnt;
    const auto etaij = getcol(eta,i,j);
    const auto g = [&] (const int k, Real& accum, const bool final) {
      assert(k != 0 or accum == 0);
      const Real d = k == 0 ? etaij(0) + detaij(0) : detaij(k);
      accum += d;
      if (final) etaij(k+1) = accum;
    };
    Dispatch<>::parallel_scan(kv.team, nlev-1, g);
  };
  int cnt = 0;
  Kokkos::parallel_reduce(ttr, f2, cnt);
  return cnt;
}

// Compute surface pressure ps = ai(0) ps0 + sum dp.
KOKKOS_FUNCTION void calc_ps (
  const KernelVariables& kv, const int nlev,
  const Real& ps0, const Real& hyai0,
  const CSelnV& dp,
  const ExecViewUnmanaged<Real[NP][NP]>& ps)
{
  assert_eln(dp, nlev);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr_snlev = Kokkos::ThreadVectorRange(kv.team, nlev);
  const CRelnV dps = elp2r(dp);
  const auto f1 = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto g = [&] (int k, Real& sum) { sum += dps(i,j,k); };
    Real sum;
    Dispatch<>::parallel_reduce(kv.team, tvr_snlev, g, sum);
    Kokkos::single(Kokkos::PerThread(kv.team),
                   [&] { ps(i,j) = hyai0*ps0 + sum; });
  };
  Kokkos::parallel_for(ttr, f1);
}

// Compute the surface pressure ps[i] at time point i corresponding to
// dp[i] = (1-alpha[i]) dp1 + alpha[i] dp2.
template <typename Snapshots>
KOKKOS_FUNCTION void calc_ps (
  const KernelVariables& kv, const int nlev,
  const Real& ps0, const Real& hyai0,
  const Snapshots& snaps,
  const ExecViewUnmanaged<Real[2][NP][NP]>& ps)
{
  const auto ie = kv.ie;
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr_snlev = Kokkos::ThreadVectorRange(kv.team, nlev);
  const auto f1 = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    for (int t = 0; t < 2; ++t) {
      const auto e = snaps.get_element(ie);
      const auto g = [&] (int k, Real& sum) {
        sum += e.get_dp_real(t,i,j,k);
      };
      Real sum;
      Dispatch<>::parallel_reduce(kv.team, tvr_snlev, g, sum);
      Kokkos::single(Kokkos::PerThread(kv.team), [&] { ps(t,i,j) = sum; });
    }
  };
  Kokkos::parallel_for(ttr, f1);
  kv.team_barrier();
  const auto f2 = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto g = [&] () {
      Real vals[2];
      for (int t = 0; t < 2; ++t) {
        const auto alpha = snaps.get_alpha(t);
        vals[t] = (hyai0*ps0 +
                   (1 - alpha)*ps(0,i,j) +
                   /**/ alpha *ps(1,i,j));
      }
      for (int t = 0; t < 2; ++t)
        ps(t,i,j) = vals[t];
    };
    Kokkos::single(Kokkos::PerThread(kv.team), g);
  };
  Kokkos::parallel_for(ttr, f2);
}

// Transform eta_dot_dpdn to eta_dot, both at interfaces, using the formula
//     eta_dot = eta_dot_dpdn/(A_eta p0 + B_eta ps)
//             = eta_dot_dpdn/(p0 + B_eta (ps - p0)).
KOKKOS_FUNCTION void calc_etadotint_from_etadotdpdnint (
  const KernelVariables& kv, const int nlev,
  const Real ps0, const CSnV& db_deta_i, const CRelV& ps,
  //  in: eta_dot_dpdn at interfaces
  // out: eta_dot at interfaces
  const SelnV& ed)
{
  assert(calc_nscal(db_deta_i.extent_int(0)) >= nlev+1);
  assert_eln(ed, nlev+1);
  const auto f = [&] (const int i, const int j, const int k) {
    ed(i,j,k) = ed(i,j,k) / (ps0 + db_deta_i(k)*(ps(i,j) - ps0));
  };
  cti::loop_ijk(calc_npack(nlev), kv, f); // final level is unused
}

// Compute eta_dot at midpoint or interface nodes at the start and end of the
// substep.
template <typename Snapshots>
KOKKOS_FUNCTION void calc_eta_dot_ref (
  const KernelVariables& kv, const SphereOperators& sphops, const Snapshots& snaps,
  const Real& ps0, const Real& hyai0, const CSNV<NUM_LEV_P>& hybi,
  const CSNV<NUM_LEV>& hydai, const CSNV<NUM_LEV>& hydbi, // delta ai, bi
  const CSNV<NUM_LEV>& hydetai, // delta etai
  const CSNV<NUM_LEV_P>& db_deta_i,
  const SelNlevp& wrk1, const SelNlevp& wrk2, const S2elNlevp& vwrk1,
  // Holds interface levels as intermediate data but is midpoint data on output,
  // with final slot unused.
  const SelNlevp eta_dot[2])
{
  using Kokkos::ALL;
  const auto ie = kv.ie;
  const int nlev = NUM_PHYSICAL_LEV;
  const SelNlev divdp(wrk1.data());
  const S2elNlev vdp(vwrk1.data());
  const ExecViewUnmanaged<Real[2][NP][NP]> ps(cti::pack2real(wrk2));
  // Calc surface pressure for use at the end.
  calc_ps(kv, nlev,
          ps0, hyai0,
          snaps,
          ps);
  kv.team_barrier();
  for (int t = 0; t < 2; ++t) {
    // Compute divdp.
    const auto e = snaps.get_element(ie);
    const auto f = [&] (const int i, const int j, const int kp) {
      for (int d = 0; d < 2; ++d)
        vdp(d,i,j,kp) = e.combine_vdp(t,d,i,j,kp);
    };
    cti::loop_ijk<cti::num_lev_pack>(kv, f);
    kv.team_barrier();
    sphops.divergence_sphere(kv, vdp, divdp);
    kv.team_barrier();
    // Compute eta_dot_dpdn at interface nodes.
    const auto& edd = eta_dot[t];
    const RelNlevp edds(cti::pack2real(edd));
    const RelNlev divdps(cti::pack2real(wrk1));
    cti::calc_eta_dot_dpdn(kv,
                           hybi,
                           divdps, edd,
                           edds);
    kv.team_barrier();
    const auto pst = Kokkos::subview(ps,t,ALL,ALL);
    calc_etadotint_from_etadotdpdnint(kv, nlev, ps0, db_deta_i, pst, edd);
    // No team_barrier: wrk1 is protected in second iteration.
  }
}

// Given the vertical and horizontal nodal velocities at time endpoints,
// evaluate the velocity estimate formula, providing the final horizontal
// velocity estimates at midpoint nodes.
KOKKOS_FUNCTION void calc_vel_horiz_formula_node_ref_mid (
  const KernelVariables& kv, const SphereOperators& sphere_ops,
  const CSNV<NUM_LEV>& hyetam, const ExecViewUnmanaged<Real[2][3][NP][NP]>& vec_sph2cart,
  // Velocities are at midpoints. Final eta_dot entry is ignored.
  const Real dtsub, const CS2elNlev vsph[2], const CSelNlevp eta_dot[2],
  const SelNlevp& wrk1, const S2elNlevp& vwrk1, const S2elNlevp& vwrk2,
  const ExecViewUnmanaged<Real****>& vnode)
{
  using Kokkos::ALL;
  const S2elNlev vfsph(vwrk1.data()), vw2(vwrk2.data());
  const SelNlev w1(wrk1.data());
  const R2elNlev vfsphs(cti::pack2real(vfsph));
  const auto& vsph1 = vsph[0];
  const auto& vsph2 = vsph[1];
  { // Horizontal terms.
    cti::ugradv_sphere(sphere_ops, kv, vec_sph2cart, vsph2, vsph1, w1, vw2, vfsph);
    for (int d = 0; d < 2; ++d) {
      const auto f = [&] (const int i, const int j, const int k) {
        vfsph(d,i,j,k) = vsph1(d,i,j,k) + vsph2(d,i,j,k) - dtsub*vfsph(d,i,j,k);
      };
      cti::loop_ijk<cti::num_lev_pack>(kv, f);
    }
  }
  kv.team_barrier();
  { // Vertical terms.
    const CRNV<NUM_PHYSICAL_LEV> etams(cti::cpack2real(hyetam));
    const CR2elNlev vsph1s(cti::cpack2real(vsph1));
    const CRelNlevp eds(cti::cpack2real(eta_dot[1]));
    for (int d = 0; d < 2; ++d) {
      const auto f = [&] (const int i, const int j, const int k) {
        Real deriv;
        if (k == 0 or k+1 == NUM_PHYSICAL_LEV) {
          const int k1 = k == 0 ? 0 : NUM_PHYSICAL_LEV-2;
          const int k2 = k == 0 ? 1 : NUM_PHYSICAL_LEV-1;
          deriv = ((vsph1s(d,i,j,k2) - vsph1s(d,i,j,k1)) /
                   (etams(k2) - etams(k1)));
        } else {
          deriv = cti::approx_derivative(
            etams(k-1), etams(k), etams(k+1),
            vsph1s(d,i,j,k-1), vsph1s(d,i,j,k), vsph1s(d,i,j,k+1));
        }
        // Interpolate eta_dot at interfaces to midpoints. Note that this is the
        // only time this is done, and it's used only in a term of the formula
        // that contains dtsub.
        const auto eta_dot = (eds(i,j,k) + eds(i,j,k+1))/2;
        vfsphs(d,i,j,k) = (vfsphs(d,i,j,k) - dtsub*eta_dot*deriv)/2;
      };
      cti::loop_ijk<cti::num_phys_lev>(kv, f);
    }
  }
  { // Transform to Cartesian.
    for (int d = 0; d < 3; ++d) {
      const auto f = [&] (const int i, const int j, const int k) {
        vnode(k,i,j,d) = (vec_sph2cart(0,d,i,j)*vfsphs(0,i,j,k) +
                          vec_sph2cart(1,d,i,j)*vfsphs(1,i,j,k));
      };
      cti::loop_ijk<cti::num_phys_lev>(kv, f);
    }
  }
}

// Given the vertical (interface) and horizontal (midpoint) nodal velocities at
// time endpoints, evaluate the velocity estimate formula, providing the final
// vertical velocity estimates at interface nodes.
KOKKOS_FUNCTION void calc_eta_dot_formula_node_ref_int (
  const KernelVariables& kv, const SphereOperators& sphere_ops,
  const CRNV<NUM_INTERFACE_LEV>& hyetai, const CSNV<NUM_LEV>& hyetam,
  // Velocities are at midpoints. Final eta_dot entry is ignored.
  const Real dtsub, const CS2elNlev vsph[2], const CSelNlevp eta_dot[2],
  const SelNlevp& wrk1, const S2elNlevp& vwrk1,
  const ExecViewUnmanaged<Real****>& vnode)
{
  const RelNlev ed1_vderiv(cti::pack2real(wrk1));
  const CRNV<NUM_INTERFACE_LEV> etai(hyetai);
  { // \dot{eta}_eta evaluated at interfaces
    const CRelNlevp ed1s(cti::cpack2real(eta_dot[0]));
    const auto f = [&] (const int i, const int j, const int km1) {
      const auto k = km1 + 1;
      ed1_vderiv(i,j,k) =
        cti::approx_derivative(
          etai(k-1), etai(k), etai(k+1),
          ed1s(i,j,k-1), ed1s(i,j,k), ed1s(i,j,k+1));
    };
    cti::loop_ijk<cti::num_phys_lev-1>(kv, f);
  }
  kv.team_barrier();
  const S2elNlev ed1_hderiv_p(vwrk1.data());
  sphere_ops.gradient_sphere(kv, eta_dot[0], ed1_hderiv_p, NUM_LEV);
  const CR2elNlev ed1_hderiv(cti::pack2real(ed1_hderiv_p));
  const CRNV<NUM_PHYSICAL_LEV> etam(cti::cpack2real(hyetam));
  {
    const CR2elNlev vsph2(cti::cpack2real(vsph[1]));
    const CRelNlevp ed1(cti::cpack2real(eta_dot[0]));
    const CRelNlevp ed2(cti::cpack2real(eta_dot[1]));
    const auto f = [&] (const int i, const int j, const int km1) {
      const auto k = km1 + 1;
      // Linearly interp horiz velocity to interfaces.
      const auto a = (etai(k) - etam(k-1)) / (etam(k) - etam(k-1));
      vnode(k,i,j,3) =
        (ed1(i,j,k) + ed2(i,j,k)
         - dtsub*(  ((1-a)*vsph2(0,i,j,k-1) + a*vsph2(0,i,j,k))*ed1_hderiv(0,i,j,k)
                  + ((1-a)*vsph2(1,i,j,k-1) + a*vsph2(1,i,j,k))*ed1_hderiv(1,i,j,k)
                  + ed2(i,j,k)*ed1_vderiv(i,j,k)))/2;
    };
    cti::loop_ijk<cti::num_phys_lev-1>(kv, f);
  }
}

} // namespace anon
} // namespace Homme

#endif
#endif
