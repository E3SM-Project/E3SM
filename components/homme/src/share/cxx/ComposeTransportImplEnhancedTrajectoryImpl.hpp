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
    if (caller)
      printf("linterp: xi out of bounds: %s %1.15e %1.15e %1.15e %1.15e\n",
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

// Compute Lagrangian level midpoints at t1 on arrival column:
//     eta_arr_mid = I[eta_ref_mid([eta(0),eta_dep_mid,eta(1)])](eta_ref_mid).
KOKKOS_FUNCTION void
eta_interp_eta (const KernelVariables& kv, const int nlev,
                const CRnV& hy_etai, const CRelnV& x, const CRnV& y,
                const RelnV& xwrk, const RnV& ywrk,
                // Use xi(i_os:), yi(i,j,i_os:).
                const int ni, const CRnV& xi, const RelnV& yi, const int i_os = 0) {
  const auto& xbdy = xwrk;
  const auto& ybdy = ywrk;
  assert(hy_etai.extent_int(0) >= nlev+1);
  assert_eln(x, nlev);
  assert(y.extent_int(0) >= nlev);
  assert_eln(xbdy, nlev+2);
  assert(ybdy.extent_int(0) >= nlev+2);
  assert(xi.extent_int(0) >= i_os + ni);
  assert_eln(yi, i_os + ni);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr_ni = Kokkos::ThreadVectorRange(kv.team, ni);
  const auto tvr_nlevp2 = Kokkos::ThreadVectorRange(kv.team, nlev+2);
  const auto f_y = [&] (const int k) {
    ybdy(k) = (k == 0      ? hy_etai(0) :
               k == nlev+1 ? hy_etai(nlev) :
               /**/          y(k-1));
  };
  Kokkos::parallel_for(Kokkos::TeamVectorRange(kv.team, nlev+2), f_y);
  kv.team_barrier();
  const auto f_x = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto g = [&] (const int k) {
      xbdy(i,j,k) = (k == 0      ? hy_etai(0) :
                     k == nlev+1 ? hy_etai(nlev) :
                     /**/          x(i,j,k-1));
    };
    Kokkos::parallel_for(tvr_nlevp2, g);
  };
  Kokkos::parallel_for(ttr, f_x);
  kv.team_barrier();
  const auto f_linterp = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    linterp(tvr_ni,
            nlev+2, getcolc(xbdy,i,j), ybdy,
            ni,     xi.data() + i_os, getcol(yi,i,j).data() + i_os,
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
deta_caas (const KernelVariables& kv, const Range& tvr_nlevp,
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
  Dispatch<>::parallel_reduce(kv.team, tvr_nlevp, g1, sums);
  const Real wneeded = sums.v[0];
  if (wneeded == 0) return;
  // Remove what is needed from the donors.
  const Real wavail = sums.v[1];
  const auto g2 = [&] (const int k) {
    deta(k) += wneeded*(w(k)/wavail);
  };
  Kokkos::parallel_for(tvr_nlevp, g2);
}

// Wrapper to above.
KOKKOS_FUNCTION void
deta_caas (const KernelVariables& kv, const int nlevp, const CRnV& deta_ref,
           const Real low, const RelnV& wrk, const RelnV& deta) {
  assert(deta_ref.extent_int(0) >= nlevp);
  assert_eln(wrk, nlevp);
  assert_eln(deta, nlevp);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr = Kokkos::ThreadVectorRange(kv.team, nlevp);
  const auto f = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    deta_caas(kv, tvr, deta_ref, low, getcol(wrk,i,j), getcol(deta,i,j));
  };
  Kokkos::parallel_for(ttr, f);
}

// Wrapper to deta_caas. On input and output, eta contains the midpoint eta
// values. On output, deta_caas has been applied, if necessary, to
// diff(eta(i,j,:)).
KOKKOS_FUNCTION void
limit_etam (const KernelVariables& kv, const int nlev, const CRnV& hy_etai,
            const CRnV& deta_ref, const Real deta_tol, const RelnV& wrk1,
            const RelnV& wrk2, const RelnV& eta) {
  assert(hy_etai.extent_int(0) >= nlev+1);
  assert(deta_ref.extent_int(0) >= nlev+1);
  const auto deta = wrk2;
  assert_eln(wrk1, nlev+1);
  assert_eln(deta, nlev+1);
  assert_eln(eta , nlev  );
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr = Kokkos::ThreadVectorRange(kv.team, nlev+1);
  // eta -> deta; limit deta if needed.
  const auto f1 = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto  etaij = getcolc( eta,i,j);
    const auto detaij = getcol(deta,i,j);
    const auto g1 = [&] (const int k, int& nbad) {
      const auto d = (k == 0    ? etaij(0) - hy_etai(0) :
                      k == nlev ? hy_etai(nlev) - etaij(nlev-1) :
                      /**/        etaij(k) - etaij(k-1));
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
    };
    deta_caas(kv, tvr, deta_ref, deta_tol, getcol(wrk1,i,j), detaij);
  };
  Kokkos::parallel_for(ttr, f1);
  kv.team_barrier();
  // deta -> eta; ignore columns where limiting wasn't needed.
  const auto f2 = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    const auto  etaij = getcol( eta,i,j);
    const auto detaij = getcol(deta,i,j);
    if (detaij(0) == -1) return;
    const auto g = [&] (const int k, Real& accum, const bool final) {
      assert(k != 0 or accum == 0);
      const Real d = k == 0 ? hy_etai(0) + detaij(0) : detaij(k);
      accum += d;
      if (final) etaij(k) = accum;
    };
    Dispatch<>::parallel_scan(kv.team, nlev, g);
  };
  Kokkos::parallel_for(ttr, f2);
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
KOKKOS_FUNCTION void calc_ps (
  const KernelVariables& kv, const int nlev,
  const Real& ps0, const Real& hyai0,
  const Real alpha[2], const CSelnV& dp1, const CSelnV& dp2,
  const ExecViewUnmanaged<Real[2][NP][NP]>& ps)
{
  assert_eln(dp1, nlev);
  assert_eln(dp2, nlev);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr_snlev = Kokkos::ThreadVectorRange(kv.team, nlev);
  const CRelnV dps[] = {elp2r(dp1), elp2r(dp2)};
  const auto f1 = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    for (int t = 0; t < 2; ++t) {
      const auto& dp = dps[t];
      const auto g = [&] (int k, Real& sum) { sum += dp(i,j,k); };
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
      for (int t = 0; t < 2; ++t)
        vals[t] = (hyai0*ps0 +
                   (1 - alpha[t])*ps(0,i,j) +
                   /**/ alpha[t] *ps(1,i,j));
      for (int t = 0; t < 2; ++t)
        ps(t,i,j) = vals[t];
    };
    Kokkos::single(Kokkos::PerThread(kv.team), g);
  };
  Kokkos::parallel_for(ttr, f2);
}

// Transform eta_dot_dpdn at interfaces to eta_dot at midpoints using the
// formula
//     eta_dot = eta_dot_dpdn/(A_eta p0 + B_eta ps).
//            a= eta_dot_dpdn diff(eta)/(diff(A) p0 + diff(B) ps).
KOKKOS_FUNCTION void calc_etadotmid_from_etadotdpdnint (
  const KernelVariables& kv, const int nlev,
  const Real& ps0, const CSnV& hydai, const CSnV& hydbi,
  const CSnV& hydetai, const CRelV& ps, const SelnV& wrk,
  //  in: eta_dot_dpdn at interfaces
  // out: eta_dot at midpoints, final slot unused
  const SelnV& ed)
{
  assert(calc_nscal(hydai.extent_int(0)) >= nlev);
  assert(calc_nscal(hydbi.extent_int(0)) >= nlev);
  assert(calc_nscal(hydetai.extent_int(0)) >= nlev);
  assert_eln(wrk, nlev+1);
  assert_eln(ed, nlev+1);
  const auto& edd_mid = wrk;
  {
    const CRelnV edd(elp2r(ed));
    const RelnV tmp(elp2r(wrk));
    const auto f = [&] (const int i, const int j, const int k) {
      tmp(i,j,k) = (edd(i,j,k) + edd(i,j,k+1))/2;
    };
    cti::loop_ijk(nlev, kv, f);
  }
  kv.team_barrier();
  {
    const auto f = [&] (const int i, const int j, const int kp) {
      ed(i,j,kp) = (edd_mid(i,j,kp)
                    * hydetai(kp)
                    / (hydai(kp)*ps0 + hydbi(kp)*ps(i,j)));
    };
    cti::loop_ijk(calc_npack(nlev), kv, f);
  }
}

// Compute eta_dot at midpoint nodes at the start and end of the substep.
KOKKOS_FUNCTION void calc_eta_dot_ref_mid (
  const KernelVariables& kv, const SphereOperators& sphere_ops,
  const Real& ps0, const Real& hyai0, const CSNV<NUM_LEV_P>& hybi,
  const CSNV<NUM_LEV>& hydai, const CSNV<NUM_LEV>& hydbi, // delta ai, bi
  const CSNV<NUM_LEV>& hydetai, // delta etai
  const Real alpha[2],
  const CS2elNlev& v1, const CSelNlev& dp1, const CS2elNlev& v2, const CSelNlev& dp2,
  const SelNlevp& wrk1, const SelNlevp& wrk2, const S2elNlevp& vwrk1,
  // Holds interface levels as intermediate data but is midpoint data on output,
  // with final slot unused.
  const SelNlevp eta_dot[2])
{
  using Kokkos::ALL;
  const int nlev = NUM_PHYSICAL_LEV;
  const SelNlev divdp(wrk1.data());
  const S2elNlev vdp(vwrk1.data());
  const ExecViewUnmanaged<Real[2][NP][NP]> ps(cti::pack2real(wrk2));
  // Calc surface pressure for use at the end.
  calc_ps(kv, nlev,
          ps0, hyai0,
          alpha, dp1, dp2,
          ps);
  kv.team_barrier();
  for (int t = 0; t < 2; ++t) {
    // Compute divdp.
    const auto f = [&] (const int i, const int j, const int kp) {
      for (int d = 0; d < 2; ++d)
        vdp(d,i,j,kp) = ((1 - alpha[t])*v1(d,i,j,kp)*dp1(i,j,kp) +
                         /**/ alpha[t] *v2(d,i,j,kp)*dp2(i,j,kp));
    };
    cti::loop_ijk<cti::num_lev_pack>(kv, f);
    kv.team_barrier();
    sphere_ops.divergence_sphere(kv, vdp, divdp);
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
    calc_etadotmid_from_etadotdpdnint(kv, nlev,
                                      ps0, hydai, hydbi, hydetai,
                                      Kokkos::subview(ps,t,ALL,ALL),
                                      wrk1,
                                      edd);
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
        vfsphs(d,i,j,k) = (vfsphs(d,i,j,k) - dtsub*eds(i,j,k)*deriv)/2;
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

// Given the vertical and horizontal nodal velocities at time endpoints,
// evaluate the velocity estimate formula, providing the final vertical velocity
// estimates at midpoint nodes.
KOKKOS_FUNCTION void calc_eta_dot_formula_node_ref_mid (
  const KernelVariables& kv, const SphereOperators& sphere_ops,
  const CRNV<NUM_INTERFACE_LEV>& hyetai, const CSNV<NUM_LEV>& hyetam,
  // Velocities are at midpoints. Final eta_dot entry is ignored.
  const Real dtsub, const CS2elNlev vsph[2], const CSelNlevp eta_dot[2],
  const SelNlevp& wrk1, const S2elNlevp& vwrk1,
  const ExecViewUnmanaged<Real****>& vnode)
{
  const SelNlev ed1_vderiv(wrk1.data());
  {
    const CRNV<NUM_PHYSICAL_LEV> etams(cti::cpack2real(hyetam));
    const CRelNlevp ed1s(cti::cpack2real(eta_dot[0]));
    const RelNlev ed1_vderiv_s(cti::pack2real(ed1_vderiv));
    const auto f = [&] (const int i, const int j, const int k) {
      Real deriv;
      if (k == 0 or k+1 == NUM_PHYSICAL_LEV) {
        deriv = cti::approx_derivative(
          k == 0 ? hyetai(0) : etams(k-1),
          etams(k),
          k+1 == NUM_PHYSICAL_LEV ? hyetai(NUM_PHYSICAL_LEV) : etams(k+1),
          k == 0 ? 0 : ed1s(i,j,k-1),
          ed1s(i,j,k),
          k+1 == NUM_PHYSICAL_LEV ? 0 : ed1s(i,j,k+1));
      } else {
        deriv = cti::approx_derivative(
          etams(k-1), etams(k), etams(k+1),
          ed1s(i,j,k-1), ed1s(i,j,k), ed1s(i,j,k+1));
      }
      ed1_vderiv_s(i,j,k) = deriv;
    };
    cti::loop_ijk<cti::num_phys_lev>(kv, f);
  }
  kv.team_barrier();
  const S2elNlev ed1_hderiv(vwrk1.data());
  sphere_ops.gradient_sphere(kv, eta_dot[0], ed1_hderiv, NUM_LEV);
  {
    const auto& vsph2 = vsph[1];
    const auto& ed1 = eta_dot[0];
    const auto& ed2 = eta_dot[1];
    const auto f = [&] (const int i, const int j, const int k) {
      const auto v = (ed1(i,j,k) + ed2(i,j,k)
                      - dtsub*(  vsph2(0,i,j,k)*ed1_hderiv(0,i,j,k)
                               + vsph2(1,i,j,k)*ed1_hderiv(1,i,j,k)
                               +   ed2(  i,j,k)*ed1_vderiv(  i,j,k)))/2;
      for (int s = 0; s < VECTOR_SIZE; ++s)
        vnode(VECTOR_SIZE*k+s, i,j,3) = v[s];
    };
    cti::loop_ijk<cti::num_lev_pack>(kv, f);
  }
}

} // namespace anon
} // namespace Homme

#endif
#endif
