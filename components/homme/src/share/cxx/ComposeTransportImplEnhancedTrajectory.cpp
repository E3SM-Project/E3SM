/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImpl.hpp"

#include "compose_hommexx.hpp"

#include <random>

namespace Homme {

// For limit_etam.
void ComposeTransportImpl::setup_enhanced_trajectory () {
  const auto etai = cmvdc(m_hvcoord.etai);
  const Real deta_ave = (etai(num_phys_lev) - etai(0)) / num_phys_lev;
  m_data.deta_tol = 10*std::numeric_limits<Real>::epsilon()*deta_ave;

  // diff(etai)
  m_data.hydetai = decltype(m_data.hydetai)("hydetai");
  const auto detai_pack = Kokkos::create_mirror_view(m_data.hydetai);
  HostViewUnmanaged<Real[NUM_PHYSICAL_LEV]> detai(pack2real(detai_pack));
  for (int k = 0; k < num_phys_lev; ++k)
    detai(k) = etai(k+1) - etai(k);
  Kokkos::deep_copy(m_data.hydetai, detai_pack);

  const auto etamp = cmvdc(m_hvcoord.etam);
  HostViewUnmanaged<Real[NUM_PHYSICAL_LEV]> etam(pack2real(etamp));
  
  // hydetam_ref.
  m_data.hydetam_ref = decltype(m_data.hydetam_ref)("hydetam_ref");
  const auto m = Kokkos::create_mirror_view(m_data.hydetam_ref);
  const int nlev = num_phys_lev;
  m(0) = etam(0) - etai(0);
  for (int k = 1; k < nlev; ++k) m(k) = etam(k) - etam(k-1);
  m(nlev) = etai(nlev) - etam(nlev-1);
  Kokkos::deep_copy(m_data.hydetam_ref, m);

  // etam
  homme::compose::set_hvcoord(etai(0), etai(num_phys_lev), etam.data());
}

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

// For sorted ascending x[0:n] and x in [x[0], x[n-1]] with hint xi_idx, return
// i such that x[i] <= xi <= x[i+1].
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

// Linear interpolation core computation.
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

// Set dep_points_all to level-midpoint arrival points.
void init_dep_points (const CTI& c, const cti::DeparturePoints& dep_pts) {
  const auto independent_time_steps = c.m_data.independent_time_steps;
  const auto& sphere_cart = c.m_geometry.m_sphere_cart;
  const CRNV<NUM_PHYSICAL_LEV> hyetam(cti::cpack2real(c.m_hvcoord.etam));
  assert(not independent_time_steps or dep_pts.extent_int(4) == 4);
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < 3; ++d)
      dep_pts(ie,lev,i,j,d) = sphere_cart(ie,i,j,d);
    if (independent_time_steps)
      dep_pts(ie,lev,i,j,3) = hyetam(lev);
  };
  c.launch_ie_physlev_ij(f);
}

void update_dep_points (
  const CTI& c, const Real dtsub, const cti::DeparturePoints& vdep,
  const cti::DeparturePoints& dep_pts)
{
  const auto independent_time_steps = c.m_data.independent_time_steps;
  const auto is_sphere = c.m_data.geometry_type == 0;
  const auto scale_factor = c.m_geometry.m_scale_factor;
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    // Update horizontal position.
    Real p[3];
    for (int d = 0; d < 3; ++d)
      p[d] = dep_pts(ie,lev,i,j,d) - dtsub*vdep(ie,lev,i,j,d)/scale_factor;
    if (is_sphere) {
      const auto norm = std::sqrt(square(p[0]) + square(p[1]) + square(p[2]));
      for (int d = 0; d < 3; ++d)
        p[d] /= norm;
    }
    for (int d = 0; d < 3; ++d)
      dep_pts(ie,lev,i,j,d) = p[d];
    if (independent_time_steps) {
      // Update vertical position.
      dep_pts(ie,lev,i,j,3) -= dtsub*vdep(ie,lev,i,j,3);
    }
  };
  c.launch_ie_physlev_ij(f);
}

/* Evaluate a formula to provide an estimate of nodal velocities that are use to
   create a 2nd-order update to the trajectory. The fundamental formula for the
   update in position p from arrival point p1 to departure point p0 is
       p0 = p1 - dt/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).
   Here we compute the velocity estimate at the nodes:
       1/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).
*/
void calc_nodal_velocities (
  const CTI& c, const Real dtsub, const Real halpha[2],
  const cti::CVSlot& v1, const cti::CDpSlot& dp1, const int idx1,
  const cti::CVSlot& v2, const cti::CDpSlot& dp2, const int idx2,
  const cti::DeparturePoints& vnode)
{
  using Kokkos::ALL;
  const auto& d = c.m_data;
  const auto& h = c.m_hvcoord;
  const auto& sphere_ops = c.m_sphere_ops;
  const auto& vec_sph2cart = c.m_geometry.m_vec_sph2cart;
  const bool independent_time_steps = d.independent_time_steps;
  const auto ps0 = h.ps0;
  const auto hyai0 = h.hybrid_ai0;
  const auto& hybi = h.hybrid_bi_packed;
  const auto& hydai = h.hybrid_ai_delta;
  const auto& hydbi = h.hybrid_bi_delta;
  const auto& hyetam = h.etam;
  const auto& hyetai = h.etai;
  const auto& hydetai = d.hydetai;
  const auto& buf1a = d.buf1o[0]; const auto& buf1b = d.buf1o[1];
  const auto& buf1c = d.buf1o[2]; const auto& buf1d = d.buf1o[3];
  const auto& buf2a = d.buf2 [0]; const auto& buf2b = d.buf2 [1];
  const auto& buf2c = d.buf2 [2]; const auto& buf2d = d.buf2 [3];
  const auto alpha0 = halpha[0], alpha1 = halpha[1];
  const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
    KernelVariables kv(team);
    const int ie = kv.ie;
    const auto  wrk1 = Homme::subview(buf1a, kv.team_idx);
    const auto  wrk2 = Homme::subview(buf1b, kv.team_idx);
    const auto vwrk1 = Homme::subview(buf2a, kv.team_idx);
    const auto vwrk2 = Homme::subview(buf2b, kv.team_idx);
    const auto v1_ie = Homme::subview(v1, ie, idx1);
    const auto v2_ie = Homme::subview(v2, ie, idx2);
    const Real alpha[] = {alpha0, alpha1};
    CSelNlevp eta_dot[] = {Homme::subview(buf1c, kv.team_idx),
                           Homme::subview(buf1d, kv.team_idx)};
    {
      SelNlevp eta_dot[] = {Homme::subview(buf1c, kv.team_idx),
                            Homme::subview(buf1d, kv.team_idx)};
      if (independent_time_steps) {
        const auto dp1_ie = Homme::subview(dp1, ie, idx1);
        const auto dp2_ie = Homme::subview(dp2, ie, idx2);
        calc_eta_dot_ref_mid(kv, sphere_ops,
                             ps0, hyai0, hybi, hydai, hydbi, hydetai,
                             alpha, v1_ie, dp1_ie, v2_ie, dp2_ie,
                             wrk1, wrk2, vwrk1,
                             eta_dot);
      } else {
        for (int t = 0; t < 2; ++t) {
          const auto& ed = eta_dot[t];
          const auto f = [&] (const int i, const int j, const int k) {
            ed(i,j,k) = 0;
          };
          cti::loop_ijk<cti::num_lev_pack>(kv, f);
        }
      }
    }
    // Collect the horizontal nodal velocities. v1,2 are on Eulerian levels. v1
    // is from time t1 < t2.
    auto* vm1 = Homme::subview(buf2c, kv.team_idx).data();
    auto* vm2 = Homme::subview(buf2d, kv.team_idx).data();
    CS2elNlev vsph[] = {CS2elNlev(vm1), CS2elNlev(vm2)};
    {
      S2elNlev vsph[] = {S2elNlev(vm1), S2elNlev(vm2)};
      for (int t = 0; t < 2; ++t) {
        const auto& v = vsph[t];
        for (int d = 0; d < 2; ++d) {
          const auto f = [&] (const int i, const int j, const int k) {
            v(d,i,j,k) = ((1 - alpha[t])*v1_ie(d,i,j,k) +
                          /**/ alpha[t] *v2_ie(d,i,j,k));
          };
          cti::loop_ijk<cti::num_lev_pack>(kv, f);
        }
      }
    }
    kv.team_barrier();
    // Given the vertical and horizontal nodal velocities at time endpoints,
    // evaluate the velocity estimate formula, providing the final horizontal
    // and vertical velocity estimates at midpoint nodes.
    const auto vnode_ie = Kokkos::subview(vnode, ie, ALL,ALL,ALL,ALL);
    const auto vec_sph2cart_ie = Homme::subview(vec_sph2cart, ie);
    calc_vel_horiz_formula_node_ref_mid(kv, sphere_ops,
                                        hyetam, vec_sph2cart_ie,
                                        dtsub, vsph, eta_dot,
                                        wrk1, vwrk1, vwrk2,
                                        vnode_ie);
    if (independent_time_steps) {
      kv.team_barrier();
      calc_eta_dot_formula_node_ref_mid(kv, sphere_ops,
                                        hyetai, hyetam,
                                        dtsub, vsph, eta_dot,
                                        wrk1, vwrk1,
                                        vnode_ie);
    }
  };
  Kokkos::parallel_for(c.m_tp_ne, f);
}

// Determine the departure points corresponding to the vertically Lagragnian
// grid's arrival midpoints, where the floating levels are those that evolve
// over the course of the full tracer time step. Also compute divdp, which holds
// the floating levels' dp values for later use in vertical remap.
void interp_departure_points_to_floating_level_midpoints (const CTI& c, const int np1) {
  using Kokkos::ALL;
  const int nlev = NUM_PHYSICAL_LEV, nlevp = nlev+1;
  const auto is_sphere = c.m_data.geometry_type == 0;
  const auto& d = c.m_data;
  const auto& h = c.m_hvcoord;
  const auto ps0 = h.ps0;
  const auto hyai0 = h.hybrid_ai0;
  const auto& hybi = h.hybrid_bi;
  const auto& hyetai = h.etai;
  const CRNV<NUM_PHYSICAL_LEV> hyetam(cti::cpack2real(h.etam));
  const auto& detam_ref = d.hydetam_ref;
  const auto deta_tol = d.deta_tol;
  const auto& dep_pts = d.dep_pts;
  const auto& dp3d = c.m_state.m_dp3d;
  const auto& buf1a = d.buf1e[0]; const auto& buf1b = d.buf1e[1];
  const auto& buf1c = d.buf1e[2]; const auto& buf1d = d.buf1e[3];
  const auto& buf2a = d.buf2[0];
  const auto f = KOKKOS_LAMBDA (const cti::MT& team) {
    KernelVariables kv(team);
    const int ie = kv.ie;
    const auto wrk1 = Homme::subview(buf1a, kv.team_idx);
    const auto wrk2 = Homme::subview(buf1b, kv.team_idx);
    const auto wrk3 = Homme::subview(buf1c, kv.team_idx);
    const auto wrk4 = Homme::subview(buf1d, kv.team_idx);
    const auto vwrk = Homme::subview(buf2a, kv.team_idx);
    // Reconstruct Lagrangian levels at t1 on arrival column:
    //     eta_arr_int = I[eta_ref_mid([0,eta_dep_mid,1])](eta_ref_int)
    const auto etam = p2rel(wrk3.data(), nlev);
    const auto f = [&] (const int i, const int j, const int k) {
      etam(i,j,k) = dep_pts(ie,k,i,j,3);
    };
    cti::loop_ijk<cti::num_phys_lev>(kv, f);
    kv.team_barrier();
    limit_etam(kv, nlev,
               hyetai, detam_ref, deta_tol,
               p2rel(wrk1.data(), nlevp), p2rel(wrk2.data(), nlevp),
               etam);
    kv.team_barrier();
    {
      // Compute eta_arr_int.
      const auto etai_arr = p2rel(wrk4.data(), nlevp);
      eta_interp_eta(kv, nlev,
                     hyetai,
                     etam, hyetam,
                     p2rel(wrk1.data(), nlev+2), RnV(cti::pack2real(wrk2), nlev+2),
                     nlevp-2, hyetai, etai_arr, 1);
      const auto f = [&] (const int i, const int j) {
        etai_arr(i,j,0) = hyetai(0);
        etai_arr(i,j,nlev) = hyetai(nlev);
      };
      c.loop_ij(kv, f);
      // Compute divdp.
      const ExecViewUnmanaged<Real[NP][NP]> ps(cti::pack2real(vwrk));
      calc_ps(kv, nlev,
              ps0, hyai0,
              Homme::subview(dp3d, ie, np1),
              ps);
      kv.team_barrier();
      eta_to_dp(kv, nlev,
                ps0, hybi, hyetai,
                ps, etai_arr,
                p2rel(wrk2.data(), nlev+1),
                RelnV(cti::pack2real(Homme::subview(c.m_derived.m_divdp, ie)),
                      NP, NP, NUM_LEV*VECTOR_SIZE));
      kv.team_barrier();
    }
    // Compute Lagrangian level midpoints at t1 on arrival column:
    //     eta_arr_mid = I[eta_ref_mid([0,eta_dep_mid,1])](eta_ref_mid)
    const auto etam_arr = p2rel(wrk4.data(), nlev);
    eta_interp_eta(kv, nlev,
                   hyetai,
                   etam, hyetam,
                   p2rel(wrk1.data(), nlev+2), RnV(cti::pack2real(wrk2), nlev+2),
                   nlev, hyetam, etam_arr);
    kv.team_barrier();
    // Compute departure horizontal points corresponding to arrival
    // Lagrangian level midpoints:
    //     p_dep_mid(eta_arr_mid) = I[p_dep_mid(eta_ref_mid)](eta_arr_mid)
    {
      const RelnV dpts_in(cti::pack2real(vwrk), NP, NP, nlev);
      const RelnV dpts_out(dpts_in.data() + NP*NP*nlev, NP, NP, nlev);
      for (int d = 0; d < 3; ++d) {
        const auto f = [&] (const int i, const int j, const int k) {
          dpts_in(i,j,k) = dep_pts(ie,k,i,j,d);
        };
        c.loop_ijk<cti::num_phys_lev>(kv, f);
        kv.team_barrier();
        eta_interp_horiz(kv, nlev,
                         hyetai,
                         hyetam, dpts_in,
                         RnV(cti::pack2real(wrk2), nlev+2), p2rel(wrk1.data(), nlev+2),
                         etam_arr, dpts_out);
        kv.team_barrier();
        const auto g = [&] (const int i, const int j, const int k) {
          dep_pts(ie,k,i,j,d) = dpts_out(i,j,k);
        };
        c.loop_ijk<cti::num_phys_lev>(kv, g);
        kv.team_barrier();
      }
      if (is_sphere) {
        // Normalize.
        const auto h = [&] (const int i, const int j, const int k) {
          Real norm = 0;
          for (int d = 0; d < 3; ++d) norm += square(dep_pts(ie,k,i,j,d));
          norm = std::sqrt(norm);
          for (int d = 0; d < 3; ++d) dep_pts(ie,k,i,j,d) /= norm;
        };
        c.loop_ijk<cti::num_phys_lev>(kv, h);
      }
    }
  };
  Kokkos::parallel_for(c.m_tp_ne, f);
}

void dss_vnode (const CTI& c, const cti::DeparturePoints& vnode) {
  const int ndim = c.m_data.independent_time_steps ? 4 : 3;
  const auto& spheremp = c.m_geometry.m_spheremp;
  const auto& rspheremp = c.m_geometry.m_rspheremp;
  const auto& vp = c.m_tracers.qtens_biharmonic;
  const ExecViewUnmanaged<Real**[NP][NP][NUM_LEV*VECTOR_SIZE]>
    v(cti::pack2real(vp), vp.extent_int(0), vp.extent_int(1));
  const auto f = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < ndim; ++d)
      v(ie,d,i,j,lev) = vnode(ie,lev,i,j,d)*spheremp(ie,i,j)*rspheremp(ie,i,j);
  };
  c.launch_ie_physlev_ij(f);
  Kokkos::fence();
  const auto be = c.m_v_dss_be[c.m_data.independent_time_steps ? 1 : 0];
  be->exchange();
  Kokkos::fence();
  const auto g = KOKKOS_LAMBDA (const int idx) {
    int ie, lev, i, j;
    cti::idx_ie_physlev_ij(idx, ie, lev, i, j);
    for (int d = 0; d < ndim; ++d)
      vnode(ie,lev,i,j,d) = v(ie,d,i,j,lev);
  };
  c.launch_ie_physlev_ij(g);
}

} // namespace anon

// When nvelocity > 2, we need to do bookkeeping to (1) accumulate intermediate
// velocity data and (2) use these data. This struct does this bookkeeping. At
// init time, it builds small integer and real arrays with indices and
// weights. During time stepping, it provides these data to the
// velocity-accumulation function and the trajectory method.
//
// Parameter short names:
//   dtf = dt_tracer_factor
//   drf = dt_remap_factor
//   nsub = semi_lagrange_trajectory_nsubstep
//   nvel = semi_lagrange_trajectory_nvelocity
struct CTI::VelocityRecord {
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

  // For n = 0:dtf, obs_slots(n,0:1) = [slot1, slot2], -1 if unused. These are
  // the slots to which velocity sample n contributes. obs_slots(0 or dtf,:) are
  // always -1.
  int obs_slots(const int n, const int k) const;

  // For n = 0:dtf, obs_wts(n,0:1) = [wt1, wt2], 0 if unused.
  Real obs_wts(const int n, const int k) const;

  // Substep end point i in 0:nsub uses velocity slots run_step(n),
  // run_step(n)-1.
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

Real& CTI::VelocityRecord::t_vel (const int i) {
  assert(_nvel > 2);
  assert(i >= 0 and i < _nvel);
  return _t_vel[i];
}

int& CTI::VelocityRecord::obs_slots (const int n, const int k) {
  assert(_nvel > 2);
  assert(n >= 0 and n <= _dtf);
  assert(k >= 0 and k <= 1);
  return _obs_slots[2*n+k];
}

Real& CTI::VelocityRecord::obs_wts (const int n, const int k) {
  assert(_nvel > 2);
  assert(n >= 0 and n <= _dtf);
  assert(k >= 0 and k <= 1);
  return _obs_wts[2*n+k];  
}

int& CTI::VelocityRecord::run_step (const int i) {
  assert(_nvel > 2);
  assert(i >= 0 and i <= _nsub);
  return _run_step[i];
}

Real CTI::VelocityRecord::t_vel (const int i) const {
  return const_cast<VelocityRecord*>(this)->t_vel(i);
}

int CTI::VelocityRecord::obs_slots (const int n, const int k) const {
  return const_cast<VelocityRecord*>(this)->obs_slots(n,k);
}

Real CTI::VelocityRecord::obs_wts (const int n, const int k) const {
  return const_cast<VelocityRecord*>(this)->obs_wts(n,k);
}

int CTI::VelocityRecord::run_step (const int i) const {
  return const_cast<VelocityRecord*>(this)->run_step(i);
}

void CTI::VelocityRecord
::init (const int dtf, const int drf_param, const int nsub, const int nvel_param) {
  const int
    drf = drf_param == 0 ? 1 : drf_param,
    navail = dtf/drf + 1,
    nvel = std::min(nvel_param == -1 ?
                    (2 + (nsub-1)/2) :  // default value
                    nvel_param,
                    std::min(nsub+1,    // can't use more than this
                             navail));  // this is the max available

  _dtf = dtf; _drf = drf; _nsub = nsub; _nvel = nvel;

  // nsub <= 1: No substepping.
  // nvel <= 2: Save velocity only at endpoints, as always occurs.
  if (nsub <= 1 or nvel <= 2) {
    _nvel = 2;
    return;
  }

  _t_vel.resize(nvel);
  _obs_slots.resize(2*(dtf+1)); _obs_wts.resize(2*(dtf+1));
  _run_step.resize(nsub+1);

  // Times at which velocity data are available.
  std::vector<int> t_avail(navail); {
    int i = 0;
    for (int n = 0; n <= dtf; ++n) {
      if (n % drf != 0) continue;
      t_avail[i] = n;
      i = i + 1;
    }
    assert(i == navail);
    assert(t_avail[navail-1] == dtf);
  }

  // Times to which we associate velocity data.
  for (int n = 0; n < nvel; ++n) {
    t_vel(n) = ((n*dtf) % (nvel-1) == 0 ?
                /**/ (n*dtf) / (nvel-1) :
                Real (n*dtf) / (nvel-1));
    assert(t_vel(n) >= 0 and t_vel(n) <= dtf);
    assert(n == 0 or t_vel(n) > t_vel(n-1));
  }

  // Build the tables mapping n in 0:dtf-1 to velocity slots to accumulate into.
  for (int n = 0; n <= dtf; ++n) {
    for (int k = 0; k < 2; ++k) {
      obs_slots(n,k) = -1;
      obs_wts(n,k) = 0;
    }
    if (n == 0 or n == dtf) continue;
    if (n % drf != 0) continue;
    const int time = n;
    int iav = -1;
    for (int i = 1; i < navail-1; ++i)
      if (time == t_avail[i]) {
        iav = i;
        break;
      }
    assert(iav > 0 and iav < navail-1);
    for (int i = 1; i < nvel-1; ++i) {
      if (t_avail[iav-1] < t_vel(i) and time > t_vel(i)) {
        obs_slots(n,0) = i;
        obs_wts(n,0) = ((t_vel(i) - t_avail[iav-1]) /
                        (t_avail[iav] - t_avail[iav-1]));
      }
      if (time <= t_vel(i) and t_avail[iav+1] > t_vel(i)) {
        obs_slots(n,1) = i;
        obs_wts(n,1) = ((t_avail[iav+1] - t_vel(i)) /
                        (t_avail[iav+1] - t_avail[iav]));
      }
    }
  }

  // Build table mapping n to interval to use. The trajectories go backward in
  // time, and this table reflects that.
  run_step(0) = nvel-1;
  run_step(nsub) = 1;
  for (int n = 1; n < nsub; ++n) {
    const auto time = Real((nsub-n)*dtf)/nsub;
    int ifnd = -1;
    for (int i = 0; i < nvel-1; ++i)
      if (t_vel(i) <= time and time <= t_vel(i+1)) {
        ifnd = i;
        break;
      }
    assert(ifnd >= 0 and ifnd < nvel-1);
    run_step(n) = ifnd + 1;
  }
}

// Public function.

void ComposeTransportImpl::calc_enhanced_trajectory (const int np1, const Real dt) {
  GPTLstart("compose_calc_enhanced_trajectory");

  const auto& dep_pts = m_data.dep_pts;
  const auto& vnode = m_data.vnode;
  const auto& vdep = m_data.vdep;

  init_dep_points(*this, dep_pts);

  const int nelemd = m_data.nelemd;
  const Real dtsub = dt / m_data.trajectory_nsubstep;
  const int nsubstep = m_data.trajectory_nsubstep;
  for (int step = 0; step < nsubstep; ++step) {
    {
      Kokkos::fence();
      GPTLstart("compose_vnode");
      const Real alpha[] = {Real(nsubstep-step-1)/nsubstep,
                            Real(nsubstep-step  )/nsubstep};
      const CVSlot v1(m_derived.m_vstar.data(), nelemd, 1);
      const CDpSlot dp1(m_derived.m_dp.data(), nelemd, 1);
      const auto& v2 = m_state.m_v;
      const auto& dp2 = m_state.m_dp3d;
      calc_nodal_velocities(*this, dtsub, alpha,
                            v1, dp1, 0, v2, dp2, np1,
                            vnode);
      Kokkos::fence();
      GPTLstop("compose_vnode");
    }

    GPTLstart("compose_v_bexchv");
    dss_vnode(*this, vnode);
    Kokkos::fence();
    GPTLstop("compose_v_bexchv");

    if (step == 0) {
      update_dep_points(*this, dtsub, vnode, dep_pts);
    } else {
      GPTLstart("compose_vdep");
      homme::compose::calc_v_departure(step, dtsub);
      Kokkos::fence();
      GPTLstop("compose_vdep");

      update_dep_points(*this, dtsub, vdep, dep_pts);
    }
  }
  Kokkos::fence();

  if (m_data.independent_time_steps) {
    GPTLstart("compose_floating_dep_pts");
    interp_departure_points_to_floating_level_midpoints(*this, np1);
    Kokkos::fence();
    GPTLstop("compose_floating_dep_pts");
  }

  GPTLstop("compose_calc_enhanced_trajectory");
}

// Testing.

namespace { // anon

Kokkos::TeamPolicy<ExecSpace>
get_test_team_policy (const int nelem, const int nlev, const int ncol=NP*NP) {
  ThreadPreferences tp;
  tp.max_threads_usable = ncol;
  tp.max_vectors_usable = nlev;
  tp.prefer_threads = true;
  tp.prefer_larger_team = true;
  return Homme::get_default_team_policy<ExecSpace>(nelem, tp);
}

struct TestData {
  std::mt19937_64 engine;
  static const Real eps;
  const ComposeTransportImpl& cti;

  TestData (const CTI& cti_, const int seed = 0)
    : cti(cti_), engine(seed == 0 ? std::random_device()() : seed)
  {}

  Real urand (const Real lo = 0, const Real hi = 1) {
    std::uniform_real_distribution<Real> urb(lo, hi);
    return urb(engine);
  }
};

// Data to deal with views of packs easily in tests.
struct ColData {
  int npack;
  ExecView<Scalar*> d;
  ExecView<Scalar*>::HostMirror h;
  ExecView<Real*>::HostMirror r;

  ColData (const std::string& name, const int nlev) {
    npack = calc_npack(nlev);
    d = decltype(d)(name, npack);
    h = Kokkos::create_mirror_view(d);
    r = decltype(r)(cti::pack2real(h), calc_nscal(npack));
  }

  void h2d () { Kokkos::deep_copy(d, h); }
};

struct ElData {
  int npack;
  ExecView<Scalar***> d;
  ExecView<Scalar***>::HostMirror h;
  ExecView<Real***>::HostMirror r;

  ElData (const std::string& name, const int nlev) {
    npack = calc_npack(nlev);
    d = decltype(d)(name, NP, NP, npack);
    h = Kokkos::create_mirror_view(d);
    r = decltype(r)(cti::pack2real(h), NP, NP, calc_nscal(npack));
  }

  void d2h () { Kokkos::deep_copy(h, d); }
  void h2d () { Kokkos::deep_copy(d, h); }
};

const Real TestData::eps = std::numeric_limits<Real>::epsilon();

int test_find_support (TestData&) {
  int ne = 0;
  const int n = 97;
  std::vector<Real> x(n);
  for (int i = 0; i < n; ++i) x[i] = -11.7 + (i*i)/n;
  const int ntest = 10000;
  for (int i = 0; i < ntest; ++i) {
    const Real xi = x[0] + (Real(i)/ntest)*(x[n-1] - x[0]);
    for (int x_idx : {0, 1, n/3, n/2, n-2, n-1}) {
      const int sup = find_support(n, x.data(), x_idx, xi);
      if (sup > n-2) ++ne;
      else if (xi < x[sup] or xi > x[sup+1]) ++ne;
    }
  }
  return ne;
}

void todev (const std::vector<Real>& h, const RnV& d) {
  assert(h.size() <= d.size());
  const auto m = Kokkos::create_mirror_view(d);
  for (size_t i = 0; i < h.size(); ++i) m(i) = h[i];
  Kokkos::deep_copy(d, m);
}

void fillcols (const int n, const Real* const h, const RelnV::HostMirror& a) {
  assert(n <= a.extent_int(2));
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (size_t k = 0; k < n; ++k)
        a(i,j,k) = h[k];  
}

void todev (const int n, const Real* const h, const RelnV& d) {
  const auto m = Kokkos::create_mirror_view(d);
  fillcols(n, h, m)  ;
  Kokkos::deep_copy(d, m);
}

void todev (const std::vector<Real>& h, const RelnV& d) {
  todev(h.size(), h.data(), d);
}

void tohost (const ExecView<const Real*>& d, std::vector<Real>& h) {
  assert(h.size() <= d.size());
  const auto m = Kokkos::create_mirror_view(d);
  Kokkos::deep_copy(m, d);
  for (size_t i = 0; i < h.size(); ++i) h[i] = m(i);
}

void run_linterp (const std::vector<Real>& x, const std::vector<Real>& y,
                  std::vector<Real>& xi, std::vector<Real>& yi) {
  const auto n = x.size(), ni = xi.size();
  assert(y.size() == n); assert(yi.size() == ni);
  // input -> device (test different sizes >= n)
  ExecView<Real*> xv("xv", n), yv("yv", n+1), xiv("xiv", ni+2), yiv("yiv", ni+3);
  todev(x, xv);
  todev(y, yv);
  todev(xi, xiv);
  // call linterp
  const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
    const auto range = Kokkos::TeamVectorRange(team, ni);
    linterp(range, n, xv, yv, ni, xiv, yiv, 0, "unittest");
  };
  Homme::ThreadPreferences tp;
  tp.max_threads_usable = 1;
  tp.max_vectors_usable = ni;
  tp.prefer_threads = false;
  tp.prefer_larger_team = true;
  const auto policy = get_test_team_policy(1, n);
  Kokkos::parallel_for(policy, f);
  Kokkos::fence();
  // output -> host
  tohost(yiv, yi);
}

void make_random_sorted (TestData& td, const int n, const Real xlo, const Real xhi,
                         std::vector<Real>& x) {
  assert(n >= 2);
  x.resize(n);
  x[0] = xlo;
  for (int i = 1; i < n-1; ++i) x[i] = td.urand(xlo, xhi);
  x[n-1] = xhi;
  std::sort(x.begin(), x.end());
}

int test_linterp (TestData& td) {
  int nerr = 0;
  { // xi == x => yi == y.
    int ne = 0;
    const int n = 30;
    std::vector<Real> x(n), y(n), xi(n), yi(n);
    make_random_sorted(td, n, -0.1, 1.2, x);
    make_random_sorted(td, n, -3, -1, y);
    for (int i = 0; i < n; ++i) xi[i] = x[i];
    run_linterp(x, y, xi, yi);
    for (int i = 0; i < n; ++i)
      if (yi[i] != y[i])
        ++ne;
    nerr += ne;
  }
  { // Reconstruct a linear function exactly.
    int ne = 0;
    const int n = 56, ni = n-3;
    const Real xlo = -1.2, xhi = 3.1;
    const auto f = [&] (const Real x) { return -0.7 + 1.3*x; };
    std::vector<Real> x(n), y(n), xi(ni), yi(ni);
    for (int trial = 0; trial < 4; ++trial) {
      make_random_sorted(td, n, xlo, xhi, x);
      make_random_sorted(td, ni,
                         xlo + (trial == 1 or trial == 3 ?  0.1 : 0),
                         xhi + (trial == 2 or trial == 3 ? -0.5 : 0),
                         xi);
      for (int i = 0; i < n; ++i) y[i] = f(x[i]);
      run_linterp(x, y, xi, yi);
      for (int i = 0; i < ni; ++i)
        if (std::abs(yi[i] - f(xi[i])) > 100*td.eps)
          ++ne;
    }
    nerr += ne;
  }
  return nerr;
}

int make_random_deta (TestData& td, const Real deta_tol, const int nlev,
                      Real* const deta) {
  int nerr = 0;
  Real sum = 0;
  for (int k = 0; k < nlev; ++k) {
    deta[k] = td.urand(0, 1) + 0.1;
    sum += deta[k];
  }
  for (int k = 0; k < nlev; ++k) {
    deta[k] /= sum;
    if (deta[k] < deta_tol) ++nerr;
  }
  return nerr;
}

int make_random_deta (TestData& td, const Real deta_tol, const RnV& deta) {
  int nerr = 0;
  const int nlev = deta.extent_int(0);
  const auto m = Kokkos::create_mirror_view(deta);
  nerr = make_random_deta(td, deta_tol, nlev, &m(0));
  Kokkos::deep_copy(deta, m);
  return nerr;  
}

int make_random_deta (TestData& td, const Real deta_tol, const RelnV& deta) {
  int nerr = 0;
  const int nlev = deta.extent_int(2);
  const auto m = Kokkos::create_mirror_view(deta);
  for (int i = 0; i < NP; ++i)
    for (int j = 0; j < NP; ++j)
      nerr += make_random_deta(td, deta_tol, nlev, &m(i,j,0));
  Kokkos::deep_copy(deta, m);
  return nerr;
}

int test_deta_caas (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {15, 128, 161}) {
    const Real deta_tol = 10*td.eps/nlev;
    const auto err = [&] (const char* lbl) {
      ++nerr;
      printf("test_deta_caa nlev %d: %s\n", nlev, lbl);
    };

    // nlev+1 deltas: deta = diff([0, etam, 1])
    ExecView<Real*> deta_ref("deta_ref", nlev+1);
    ExecView<Real***> deta("deta",NP,NP,nlev+1), wrk("wrk",NP,NP,nlev+1);
    nerr += make_random_deta(td, deta_tol, deta_ref);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run = [&] (const RelnV& deta) {
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        deta_caas(kv, nlev+1, deta_ref, deta_tol, wrk, deta);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
    };

    { // Test that if all is OK, the input is not altered.
      nerr += make_random_deta(td, deta_tol, deta);
      ExecView<Real***>::HostMirror copy("copy",NP,NP,nlev+1);
      Kokkos::deep_copy(copy, deta);
      run(deta);
      const auto m = cti::cmvdc(deta);
      bool diff = false;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k <= nlev; ++k)
            if (m(i,j,k) != copy(i,j,k))
              diff = true;
      if (diff) err("input not altered");
    }

    { // Modify one etam and test that only adjacent intervals change beyond eps.
      // nlev midpoints
      ExecView<Real*> etam_ref("etam_ref",nlev);
      const auto her = Kokkos::create_mirror_view(etam_ref);
      const auto hder = cti::cmvdc(deta_ref);
      {
        her(0) = hder(0);
        for (int k = 1; k < nlev; ++k)
          her(k) = her(k-1) + hder(k);
        Kokkos::deep_copy(etam_ref, her);
      }
      std::vector<Real> etam(nlev);
      const auto hde = Kokkos::create_mirror_view(deta);
      const auto get_idx = [&] (const int i, const int j) {
        const int idx = static_cast<int>(0.15*nlev);
        return std::max(1, std::min(nlev-2, idx+NP*i+j));
      };
      for (int trial = 0; trial < 2; ++trial) {
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j) {
            for (int k = 0; k < nlev; ++k) etam[k] = her(k);
            // Perturb one level.
            const int idx = get_idx(i,j);
            etam[idx] += trial == 0 ? 1.1 : -13.1;
            hde(i,j,0) = etam[0];
            for (int k = 1; k < nlev; ++k) hde(i,j,k) = etam[k] - etam[k-1];
            hde(i,j,nlev) = 1 - etam[nlev-1];
            // Make sure we have a meaningful test.
            Real minval = 1;
            for (int k = 0; k <= nlev; ++k) minval = std::min(minval, hde(i,j,k));
            if (minval >= deta_tol) err("meaningful test");
          }
        Kokkos::deep_copy(deta, hde);
        run(deta);
        Kokkos::deep_copy(hde, deta);
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j) {
            const int idx = get_idx(i,j);
            // Min val should be deta_tol.
            Real minval = 1;
            for (int k = 0; k <= nlev; ++k) minval = std::min(minval, hde(i,j,k));
            if (minval != deta_tol) err("min val");
            // Sum of levels should be 1.
            Real sum = 0;
            for (int k = 0; k <= nlev; ++k) sum += hde(i,j,k);
            if (std::abs(sum - 1) > tol) err("sum 1");
            // Only two deltas should be affected.
            Real maxdiff = 0;
            for (int k = 0; k <= nlev; ++k) {
              const auto diff = std::abs(hde(i,j,k) - hder(k));
              if (k == idx or k == idx+1) {
                if (diff <= deta_tol) err("2 deltas a");
              } else {
                maxdiff = std::max(maxdiff, diff);
              }
            }
            if (maxdiff > tol) err("2 deltas b");
          }
      }
    }

    { // Test generally (and highly) perturbed levels.
      const auto hde = Kokkos::create_mirror_view(deta);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          Real sum = 0;
          for (int k = 0; k <= nlev; ++k) {
            hde(i,j,k) = td.urand(-0.5, 0.5);
            sum += hde(i,j,k);
          }
          // Make the column sum to 0.2 for safety in the next step.
          const Real colsum = 0.2;
          for (int k = 0; k <= nlev; ++k) hde(i,j,k) += (colsum - sum)/(nlev+1);
          for (int k = 0; k <= nlev; ++k) hde(i,j,k) /= colsum;
          sum = 0;
          for (int k = 0; k <= nlev; ++k) sum += hde(i,j,k);
          if (std::abs(sum - 1) > 10*tol) err("general sum 1");
        }
      Kokkos::deep_copy(deta, hde);
      run(deta);
      Kokkos::deep_copy(hde, deta);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          Real sum = 0, minval = 1;
          for (int k = 0; k <= nlev; ++k) sum += hde(i,j,k);
          for (int k = 0; k <= nlev; ++k) minval = std::min(minval, hde(i,j,k));
          if (std::abs(sum - 1) > 1e3*td.eps) ++nerr;
          if (minval != deta_tol) err("general minval");
        }
    }
  }
  
  return nerr;
}

struct HybridLevels {
  Real ps0, a_eta, b_eta;
  std::vector<Real> ai, dai, bi, dbi, am, bm, etai, detai, etam, detam;
};

// Follow DCMIP2012 3D tracer transport specification for a, b, eta.
void fill (HybridLevels& h, const int n) {
  h.ai.resize(n+1); h.bi.resize(n+1);
  h.am.resize(n  ); h.bm.resize(n  );
  h.etai.resize(n+1); h.etam.resize(n);

  const auto Rd = PhysicalConstants::Rgas;
  const auto T0 = 300; // K
  const auto p0 = PhysicalConstants::p0;
  const auto g = PhysicalConstants::g;
  const Real ztop = 12e3; // m

  h.ps0 = p0;

  const auto calc_pressure = [&] (const Real z) {
    return p0*std::exp(-g*z/(Rd*T0));
  };

  const Real eta_top = calc_pressure(ztop)/p0;
  assert(eta_top > 0);
  for (int i = 0; i <= n; ++i) {
    const auto z = (Real(n - i)/n)*ztop;
    h.etai[i] = calc_pressure(z)/p0;
    h.bi[i] = i == 0 ? 0 : (h.etai[i] - eta_top)/(1 - eta_top);
    h.ai[i] = h.etai[i] - h.bi[i];
    assert(i == 0 or h.etai[i] > h.etai[i-1]);
  }
  assert(h.bi  [0] == 0); // Real(n - i)/n is exactly 1, so exact = holds
  assert(h.bi  [n] == 1); // exp(0) is exactly 0, so exact = holds
  assert(h.etai[n] == 1); // same
  // b = (eta - eta_top)/(1 - eta_top) => b_eta = 1/(1 - eta_top)
  // a = eta - b => a_eta = 1 - b_eta = -eta_top/(1 - eta_top)
  // p_eta = a_eta p0 + b_eta ps
  h.b_eta = 1/(1 - eta_top);
  h.a_eta = 1 - h.b_eta;

  const auto tomid = [&] (const std::vector<Real>& in, std::vector<Real>& mi) {
    for (int i = 0; i < n; ++i) mi[i] = (in[i] + in[i+1])/2;
  };
  tomid(h.ai, h.am);
  tomid(h.bi, h.bm);
  tomid(h.etai, h.etam);

  const auto diff = [&] (const std::vector<Real>& ai, std::vector<Real>& dai) {
    dai.resize(n);
    for (int i = 0; i < n; ++i) dai[i] = ai[i+1] - ai[i];
  };
  diff(h.ai, h.dai);
  diff(h.bi, h.dbi);
  diff(h.etai, h.detai);

  h.detam.resize(n+1);
  h.detam[0] = h.etam[0] - h.etai[0];
  for (int i = 1; i < n; ++i) h.detam[i] = h.etam[i] - h.etam[i-1];
  h.detam[n] = h.etai[n] - h.etam[n-1];
}

int test_limit_etam (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {143, 128, 81}) {
    const Real deta_tol = 1e5*td.eps/nlev;

    ExecView<Real*> hy_etai("hy_etai",nlev+1), detam("detam",nlev+1);
    ExecView<Real***> wrk1("wrk1",NP,NP,nlev+1), wrk2("wrk2",NP,NP,nlev+1);
    ExecView<Real***> etam("etam",NP,NP,nlev);

    HybridLevels h;
    fill(h, nlev);
    todev(h.etai, hy_etai);
    todev(h.detam, detam);

    const auto he = Kokkos::create_mirror_view(etam);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run = [&] () {
      Kokkos::deep_copy(etam, he);
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        limit_etam(kv, nlev, hy_etai, detam, deta_tol, wrk1, wrk2, etam);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      Kokkos::deep_copy(he, etam);
    };

    fillcols(h.etam.size(), h.etam.data(), he);
    // Col 0 should be untouched. Cols 1 and 2 should have very specific changes.
    const int col1_idx = static_cast<int>(0.25*nlev);
    he(0,1,col1_idx) += 0.3;
    const int col2_idx = static_cast<int>(0.8*nlev);
    he(0,2,col2_idx) -= 5.3;
    // The rest of the columns get wild changes.
    for (int idx = 3; idx < NP*NP; ++idx) {
      const int i = idx / NP, j = idx % NP;
      for (int k = 0; k < nlev; ++k)
        he(i,j,k) += td.urand(-1, 1)*(h.etai[k+1] - h.etai[k]);
    }
    run();
    bool ok = true;
    for (int k = 0; k < nlev; ++k)
      if (he(0,0,k) != h.etam[k]) ok = false;
    for (int k = 0; k < nlev; ++k) {
      if (k == col1_idx) continue;
      if (std::abs(he(0,1,k) - h.etam[k]) > tol) ok = false;
    }
    for (int k = 0; k < nlev; ++k) {
      if (k == col2_idx) continue;
      if (std::abs(he(0,2,k) - h.etam[k]) > tol) ok = false;
    }
    Real mingap = 1;
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j) {
        mingap = std::min(mingap, he(i,j,0) - h.etai[0]);
        for (int k = 1; k < nlev; ++k)
          mingap = std::min(mingap, he(i,j,k) - he(i,j,k-1));
        mingap = std::min(mingap, h.etai[nlev] - he(i,j,nlev-1));
      }
    // Test minimum level delta, with room for numerical error.
    if (mingap < 0.8*deta_tol) ok = false;
    if (not ok) ++nerr;
  }
  
  return nerr;
}

int test_eta_interp (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {15, 128, 161}) {
    HybridLevels h;
    fill(h, nlev);

    ExecView<Real*> hy_etai("hy_etai",nlev+1);
    ExecView<Real***> x("x",NP,NP,nlev), y("y",NP,NP,nlev);
    ExecView<Real***> xi("xi",NP,NP,nlev+1), yi("yi",NP,NP,nlev+1);
    ExecView<Real***> xwrk("xwrk",NP,NP,nlev+2), ywrk("ywrk",NP,NP,nlev+2);

    todev(h.etai, hy_etai);

    const auto xh  = Kokkos::create_mirror_view(x );
    const auto yh  = Kokkos::create_mirror_view(y );
    const auto xih = Kokkos::create_mirror_view(xi);
    const auto yih = Kokkos::create_mirror_view(yi);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run_eta = [&] (const int ni) {
      Kokkos::deep_copy(x, xh); Kokkos::deep_copy(y, yh);
      Kokkos::deep_copy(xi, xih);
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        eta_interp_eta(kv, nlev, hy_etai,
                       x, getcolc(y,0,0),
                       xwrk, getcol(ywrk,0,0),
                       ni, getcolc(xi,0,0), yi);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      Kokkos::deep_copy(yih, yi);
    };
    const auto run_horiz = [&] () {
      Kokkos::deep_copy(x, xh); Kokkos::deep_copy(y, yh);
      Kokkos::deep_copy(xi, xih);
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        eta_interp_horiz(kv, nlev, hy_etai,
                         getcolc(x,0,0), y,
                         getcol(xwrk,0,0), ywrk,
                         xi, yi);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      Kokkos::deep_copy(yih, yi);
    };

    std::vector<Real> v;
    const Real d = 1e-6, vlo = h.etai[0]+d, vhi = h.etai[nlev]-d;

    for (const int ni : {int(0.7*nlev), nlev-1, nlev, nlev+1}) {
      make_random_sorted(td, nlev, vlo, vhi, v);
      fillcols(nlev, v.data(), xh);
      fillcols(nlev, v.data(), yh);
      make_random_sorted(td, ni, vlo, vhi, v);
      fillcols(ni, v.data(), xih);
      run_eta(ni);
      bool ok = true;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < ni; ++k)
            if (std::abs(yih(i,j,k) - xih(i,j,k)) > tol)
              ok = false;
      if (not ok) ++nerr;
    }

    { // Test exact interp of line in the interior, const interp near the bdys.
      make_random_sorted(td, nlev, vlo+0.05, vhi-0.1, v);
      fillcols(nlev, v.data(), xh);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          for (int k = 0; k < nlev; ++k)
            yh(i,j,k) = i*xh(0,0,k) - j;
          make_random_sorted(td, nlev, vlo, vhi, v);
          for (int k = 0; k < nlev; ++k)
            xih(i,j,k) = v[k];
        }
      run_horiz();
      bool ok = true;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k) {
            if (xih(i,j,k) < xh(0,0,0)) {
              if (std::abs(yih(i,j,k) - yih(i,j,0)) > tol)
                ok = false;
            } else if (xih(i,j,k) > xh(0,0,nlev-1)) {
              if (std::abs(yih(i,j,k) - yih(i,j,nlev-1)) > tol)
                ok = false;            
            } else {
              if (std::abs(yih(i,j,k) - (i*xih(i,j,k) - j)) > tol)
                ok = false;
            }
          }
      if (not ok) ++nerr;
    }
  }
  
  return nerr;
}

int test_eta_to_dp (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {143, 128, 81}) {
    const auto err = [&] (const char* lbl) {
      ++nerr;
      printf("test_eta_to_dp nlev %d: %s\n", nlev, lbl);
    };

    HybridLevels h;
    fill(h, nlev);

    ExecView<Real*> hy_bi("hy_bi",nlev+1), hy_etai("hy_etai",nlev+1);
    ExecView<Real***> etai("etai",NP,NP,nlev+1), wrk("wrk",NP,NP,nlev+1);
    ExecView<Real***> dp("dp",NP,NP,nlev);
    ExecView<Real[NP][NP]> ps("ps");
    const Real hy_ps0 = h.ps0;

    todev(h.bi, hy_bi);
    todev(h.etai, hy_etai);

    const auto psm = Kokkos::create_mirror_view(ps);
    HostView<Real***> dp1("dp1",NP,NP,nlev);
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j)
        psm(i,j) = (1 + 0.1*td.urand(-1, 1))*h.ps0;
    Kokkos::deep_copy(ps, psm);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run = [&] () {
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        eta_to_dp(kv, nlev, hy_ps0, hy_bi, hy_etai, ps, etai, wrk, dp);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
    };

    { // Test that for etai_ref we get the same as the usual formula.
      todev(h.etai, etai);
      HostView<Real***> dp1("dp1",NP,NP,nlev);
      Real dp1_max = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k) {
            dp1(i,j,k) = ((h.ai[k+1] - h.ai[k])*h.ps0 +
                          (h.bi[k+1] - h.bi[k])*psm(i,j));
            dp1_max = std::max(dp1_max, std::abs(dp1(i,j,k)));
          }
      run();
      const auto dph = cti::cmvdc(dp);
      Real err_max = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k)
            err_max = std::max(err_max, std::abs(dph(i,j,k) - dp1(i,j,k)));
      if (err_max > tol*dp1_max) err("t1");
    }

    { // Test that sum(dp) = ps for random input etai.
      std::vector<Real> etai_r;
      make_random_sorted(td, nlev+1, h.etai[0], h.etai[nlev], etai_r);
      todev(etai_r, etai);
      run();
      const auto dph1 = cti::cmvdc(dp);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          Real ps = h.ai[0]*h.ps0;
          for (int k = 0; k < nlev; ++k)
            ps += dph1(i,j,k);
          if (std::abs(ps - psm(i,j)) > tol*psm(i,j)) err("t2");
        }    
      // Test that values on input don't affect solution.
      Kokkos::deep_copy(wrk, 0);
      Kokkos::deep_copy(dp, 0);
      run();
      const auto dph2 = cti::cmvdc(dp);
      bool alleq = true;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k)
            if (dph2(i,j,k) != dph1(i,j,k))
              alleq = false;
      if (not alleq) err("t3");
    }
  }

  return nerr;
}

int test_calc_ps (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {15, 128, 161}) {
    HybridLevels h;
    fill(h, nlev);
    const auto ps0 = h.ps0, hyai0 = h.ai[0];

    ElData dp1("dp1", nlev), dp2("dp2", nlev);
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j)
        for (int k = 0; k < nlev; ++k) {
          dp1.r(i,j,k) = td.urand(0, 1000);
          dp2.r(i,j,k) = td.urand(0, 1000);
        }
    dp1.h2d();
    dp2.h2d();

    const Real alpha[] = {td.urand(0,1), td.urand(0,1)};

    ExecView<Real[NP][NP]> ps("ps");
    ExecView<Real[2][NP][NP]> ps2("ps2");
    const auto policy = get_test_team_policy(1, nlev);
    const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
      KernelVariables kv(team);
      calc_ps(kv, nlev, ps0, hyai0, alpha, dp1.d, dp2.d, ps2);
      calc_ps(kv, nlev, ps0, hyai0, dp1.d, ps);
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();

    const auto ps_h = cti::cmvdc(ps);
    const auto ps2_h = cti::cmvdc(ps2);
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j) {
        {
          Real ps = h.ai[0]*h.ps0;
          for (int k = 0; k < nlev; ++k)
            ps += dp1.r(i,j,k);
          if (std::abs(ps_h(i,j) - ps) > tol*ps) ++nerr;
        }
        for (int t = 0; t < 2; ++t) {
          Real ps = h.ai[0]*h.ps0;
          for (int k = 0; k < nlev; ++k)
            ps += (1 - alpha[t])*dp1.r(i,j,k) + alpha[t]*dp2.r(i,j,k);
          if (std::abs(ps2_h(t,i,j) - ps) > tol*ps) ++nerr;
        }
      }
  }

  return nerr;
}

int test_calc_etadotmid_from_etadotdpdnint (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {143, 128, 81}) {
    HybridLevels h;
    fill(h, nlev);

    // Test function:
    //     eta_dot_dpdn(eta) = c eta + d.
    // Then
    //     eta_dot = eta_dot_dpdn(eta)/dpdn(eta)
    //             = (c eta + d)/(a_eta p0 + b_eta ps).
    // Since a_eta, b_eta are constants independent of eta in this test, eta_dot
    // is then also a linear function of eta. Thus, we can test for exact
    // agreement with the true solution.

    ColData hydai("hydai",nlev), hydbi("hydbi",nlev), hydetai("hydetai",nlev);
    ElData wrk("wrk",nlev+1), ed("ed",nlev+1);
    ExecView<Real[NP][NP]> ps("ps");
    const Real ps0 = h.ps0;

    const auto ps_m = Kokkos::create_mirror_view(ps);
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j) {
        ps_m(i,j) = td.urand(0.5, 1.2)*ps0;
        for (int k = 0; k < nlev; ++k) {
          hydai.r[k] = h.dai[k];
          hydbi.r[k] = h.dbi[k];
          hydetai.r[k] = h.detai[k];
        }
        for (int k = 0; k <= nlev; ++k)
          ed.r(i,j,k) = (i-j)*h.etai[k] + 0.3;
      }
    Kokkos::deep_copy(ps, ps_m);
    hydai.h2d(); hydbi.h2d(); hydetai.h2d();
    ed.h2d();

    const auto policy = get_test_team_policy(1, nlev);
    const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
      KernelVariables kv(team);
      calc_etadotmid_from_etadotdpdnint(
        kv, nlev, ps0, hydai.d, hydbi.d, hydetai.d, ps, wrk.d, ed.d);
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    ed.d2h();

    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j) {
        const auto den = h.a_eta*h.ps0 + h.b_eta*ps_m(i,j);
        for (int k = 0; k < nlev; ++k) {
          const auto ed_true = ((i-j)*h.etam[k] + 0.3)/den;
          if (std::abs(ed.r(i,j,k) - ed_true) > tol*(10/den)) ++nerr;
        }
      }
  }

  return nerr;
}

int test_calc_eta_dot_ref_mid (TestData& td) {
  int nerr = 0;

  // calc_eta_dot_ref_mid calls several routines that are all tested
  // mathematically. calc_eta_dot_ref_mid itself is too complicated to test
  // mathematically. But we can still test it for s/w properties like
  // determinism.

  //todo

  return nerr;
}

int test_interp_departure_points_to_floating_level_midpoints (TestData& td) {
  int nerr = 0;
  //todo test case of ed = 0
  return nerr;
}

static int test1_init_velocity_record (
  const int dtf, const int drf, const int nsub, const int nvel)
{
  const auto eps = std::numeric_limits<Real>::epsilon();
  int e = 0;

  if (dtf % drf != 0) {
    printf("Testing erro: dtf %% drf == 0 is required: %d %d\n", dtf, drf);
    ++e; 
  }

  const CTI::VelocityRecord v(dtf, drf, nsub, nvel);
  if (v.dtf() != dtf) ++e;
  if (v.nsub() != nsub) ++e;

  // Check that t_vel is monotonically increasing.
  for (int n = 1; n < v.nvel(); ++n)
    if (v.t_vel(n) <= v.t_vel(n-1))
      ++e;

  // Check that obs_slots does not reference end points. This should not happen
  // b/c nvel <= navail and observations are uniformly spaced.
  for (int n = 0; n < dtf; ++n)
    for (int i = 0; i < 2; ++i)
      if (v.obs_slots(n,i) == 0 or v.obs_slots(n,i) == dtf)
        ++e;

  // Check that weights sum to 1.
  std::vector<Real> ys(dtf);
  for (int n = 0; n < dtf; ++n)
    for (int i = 0; i < 2; ++i)
      if (v.obs_slots(n,i) >= 0)
        ys[v.obs_slots(n,i)] += v.obs_wts(n,i);
  for (int i = 1; i < v.nvel()-1; ++i)
    if (std::abs(ys[i] - 1) > 1e3*eps)
      ++e;

  // Test for exact interp of an affine function.
  const auto tfn = [] (const Real x) { return 7.1*x - 11.5; };
  //   Observe data forward in time.
  Real endslots[2];
  endslots[0] = tfn(0);
  endslots[1] = tfn(dtf);
  ys[0] = -1000; // unused
  for (int i = 1; i < dtf; ++i) ys[i] = 0;
  for (int n = 1; n < dtf; ++n) {
    if (n % drf != 0) continue;
    const Real y = tfn(n);
    for (int i = 0; i < 2; ++i) {
      if (v.obs_slots(n,i) < 0) continue;
      ys[v.obs_slots(n,i)] += v.obs_wts(n,i)*y;
    }
  }
  //   Use the data backward in time.
  for (int n = 0; n < nsub; ++n) {
    // Each segment orders the data forward in time. Thus, data are always
    // ordered forward in time but used backward.
    Real xsup[2], ysup[2];
    for (int i = 0; i < 2; ++i) {
      int k = nsub - (n+1) + i;
      xsup[i] = (k*v.t_vel(v.nvel()-1))/nsub;
      k = v.run_step(nsub-i);
      const Real
        y0 = k == 1 ? endslots[0] : ys[k-1],
        y1 = k == v.nvel()-1 ? endslots[2] : ys[k];
      ysup[i] = (((v.t_vel(k) - xsup[i])*y0 + (xsup[i] - v.t_vel(k-1))*y1) /
                 (v.t_vel(k) - v.t_vel(k-1)));
    }
    for (int i = 0; i <= 10; ++i) {
      const Real
        a = Real(i)/10,
        x = (1-a)*xsup[0] + a*xsup[1],
        y = (1-a)*ysup[0] + a*ysup[1];
      if (std::abs(y - tfn(x)) > 1e3*eps) {
        printf("n %d i %2d x %7.3f y %7.3f t %7.3f\n", n, i, x, y, tfn(x));
        ++e;
      }
    }
  }
  
  if (e) {
    printf("ERROR e %d\n", e);
    printf("dtf %d drf %d nsub %d nvel %d v.nvel %d\n",
           dtf, drf, nsub, nvel, v.nvel());
    printf("  t_vel:");
    for (int i = 0; i < v.nvel(); ++i) printf(" %1.3f", v.t_vel(i));
    printf("\n  obs:\n");
    for (int n = 0; n <= dtf; ++n)
      printf("    %2d %2d %2d %1.3f %1.3f\n",
             n, v.obs_slots(n,0), v.obs_slots(n,1), v.obs_wts(n,0),
             v.obs_wts(n,1));
    printf("  run_step:\n");
    for (int n = 0; n <= nsub; ++n) printf("    %2d %2d\n", n, v.run_step(n));
  }

  return e;
}

int test_init_velocity_record (TestData& td) {
  int dtf, drf, nsub, nvel, nerr;

  nerr = 0;

  const auto f = [&] () {
    const int e = test1_init_velocity_record(dtf, drf, nsub, nvel);
    if (e > 0) ++nerr;
  };

  nerr = 0;
  dtf = 6;
  drf = 2;
  nsub = 3;
  nvel = 4;
  f();
  nvel = 3;
  f();
  drf = 3;
  nvel = 6;
  f();
  drf = 1;
  nsub = 5;
  f();
  dtf = 12;
  drf = 2;
  nsub = 3;
  nvel = -1;
  f();
  nsub = 5;
  nvel = 5;
  f();
  dtf = 27;
  drf = 3;
  nsub = 51;
  nvel = 99;
  f();

  return nerr;
}

} // namespace anon

#define comunittest(f) do {                     \
    ne = f(td);                                 \
    if (ne) printf(#f " ne %d\n", ne);          \
    nerr += ne;                                 \
  } while (0)

int ComposeTransportImpl::run_enhanced_trajectory_unit_tests () {
  int nerr = 0, ne;
  TestData td(*this);
  comunittest(test_find_support);
  comunittest(test_linterp);
  comunittest(test_eta_interp);
  comunittest(test_eta_to_dp);
  comunittest(test_deta_caas);
  comunittest(test_limit_etam);
  comunittest(test_calc_ps);
  comunittest(test_calc_etadotmid_from_etadotdpdnint);
  comunittest(test_calc_eta_dot_ref_mid);
  comunittest(test_interp_departure_points_to_floating_level_midpoints);
  comunittest(test_init_velocity_record);
  return nerr;
}

#undef comunittest

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
