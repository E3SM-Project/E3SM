/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_GLLFVREMAP_IMPL_HPP
#define HOMMEXX_GLLFVREMAP_IMPL_HPP

#include "GllFvRemap.hpp"

#include "Context.hpp"
#include "Elements.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsDerivedState.hpp"
#include "FunctorsBuffersManager.hpp"
#include "ErrorDefs.hpp"
#include "HommexxEnums.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "SphereOperators.hpp"
#include "Tracers.hpp"
#include "TimeLevel.hpp"
#include "profiling.hpp"
#include "mpi/BoundaryExchange.hpp"
#include "mpi/MpiBuffersManager.hpp"
#include "mpi/Connectivity.hpp"

#include <cassert>

namespace Homme {

struct GllFvRemapImpl {
  enum : int { np = NP };
  enum : int { packn = VECTOR_SIZE };
  enum : int { np2 = NP*NP };
  enum : int { num_lev_pack = NUM_LEV };
  enum : int { max_num_lev_pack = NUM_LEV_P };
  enum : int { num_lev_aligned = num_lev_pack*packn };
  enum : int { num_levp_aligned = max_num_lev_pack*packn };
  enum : int { num_phys_lev = NUM_PHYSICAL_LEV };

  typedef GllFvRemap::Phys1T Phys1T;
  typedef GllFvRemap::Phys2T Phys2T;
  typedef GllFvRemap::Phys3T Phys3T;
  typedef GllFvRemap::CPhys2T CPhys2T;
  typedef GllFvRemap::CPhys3T CPhys3T;
  typedef ExecViewUnmanaged<Scalar***>  VPhys2T;
  typedef ExecViewUnmanaged<Scalar****> VPhys3T;
  typedef VPhys2T::const_type CVPhys2T;
  typedef VPhys3T::const_type CVPhys3T;

  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MT = typename TeamPolicy::member_type;
  template <typename Data> using EVU = ExecViewUnmanaged<Data>;

  using Buf1 = ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV_P]>;
  using Buf2 = ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV_P]>;

  struct Data {
    int nelemd, qsize, nf2, n_dss_fld;
    bool use_moisture, theta_hydrostatic_mode;

    static constexpr int nbuf1 = 2, nbuf2 = 1;
    Buf1 buf1[nbuf1];
    Buf2 buf2[nbuf2];

    Real w_ff;
    ExecView<Real**>
      fv_metdet,   // (nelemd,nf2)
      g2f_remapd,  // (nf2,np2)
      f2g_remapd;  // (np2,nf2)
    ExecView<Real****>
      D, Dinv,     // (nelemd,np2,2,2)
      D_f, Dinv_f; // (nelemd,nf2,2,2)

    Data ()
      : nelemd(-1), qsize(-1), nf2(-1)
    {}
  };

  HybridVCoord m_hvcoord;
  Elements m_elements;
  ElementsState m_state;
  ElementsDerivedState m_derived;
  ElementsForcing m_forcing;
  ElementsGeometry m_geometry;
  Tracers m_tracers;
  Data m_data;

  TeamPolicy m_tp_ne, m_tp_ne_qsize, m_tp_ne_dss;
  TeamUtils<ExecSpace> m_tu_ne, m_tu_ne_qsize, m_tu_ne_dss;

  std::shared_ptr<BoundaryExchange> m_extrema_be, m_dss_be;

  GllFvRemapImpl();
  void setup();

  KOKKOS_INLINE_FUNCTION
  size_t shmem_size (const int team_size) const {
    return KernelVariables::shmem_size(team_size);
  }

  void reset(const SimulationParams& params);
  int requested_buffer_size() const;
  void init_buffers(const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges();

  void init_data(const int nf, const int nf_max, const bool theta_hydrostatic_mode,
                 const Real* fv_metdet_r, const Real* g2f_remapd_r,
                 const Real* f2g_remapd_r, const Real* D_f_r, const Real* Dinv_f_r);

  void run_dyn_to_fv_phys(const int time_idx, const Phys1T& ps, const Phys1T& phis,
                          const Phys2T& T, const Phys2T& omega, const Phys3T& uv,
                          const Phys3T& q, const Phys2T* dp);
  void run_fv_phys_to_dyn(const int time_idx, const CPhys2T& T, const CPhys3T& uv,
                          const CPhys3T& q);
  void run_fv_phys_to_dyn_dss();

  void remap_tracer_dyn_to_fv_phys(const int time_idx, const int nq,
                                   const CPhys3T& q_dyn, const Phys3T& q_fv);

  /* Compute pressure level increments on the FV grid given ps on the FV grid.
     Directly projecting dp_gll to dp_fv disagrees numerically with the loop in
     this subroutine. This loop is essentially how CAM computes pdel in
     derived_phys in dp_coupling.F90, so we must use it, too.
   */
  template <typename PST, typename DPT>
  static KOKKOS_FUNCTION void
  calc_dp_fv (const MT& team,
              const HybridVCoord& hvcoord, const int ncol, const int nlev,
              const PST& ps, const DPT& dp_fv) {
    assert(ps.extent_int(0) >= ncol);
    assert(dp_fv.extent_int(0) >= ncol && dp_fv.extent_int(1) >= nlev);
    using Kokkos::parallel_for;
    const auto ttr = Kokkos::TeamThreadRange(team, ncol);
    const auto tvr = Kokkos::ThreadVectorRange(team, nlev);
    loop_ik(ttr, tvr, [&] (int i, int k) {
      dp_fv(i,k) = (hvcoord.hybrid_ai_delta(k)*hvcoord.ps0 +
                    hvcoord.hybrid_bi_delta(k)*ps(i));
    });
  }

  /* Compute (1-based indexing)
         y(1:m,k) = (A (d1 x(1:n,k)))/(s2*d2), k = 1:nlev
     Sizes are min; a dim can have larger size.
         A m by n, d1 n, d2 m
         x n by nlev, w n by nlev, y m by nlev
     Permitted aliases:
         w = x
   */
  template <typename AT, typename D1T, typename D2T, typename XT, typename WT, typename YT>
  static KOKKOS_FUNCTION void
  remapd (const MT& team, const int m, const int n, const int nlev, // range of x,y fastest dim
          const AT& A, const D1T& d1, const Real s2, const D2T& d2,
          const XT& x, const YT& w, const WT& y) {
    assert(A.extent_int(0) >= m && A.extent_int(1) >= n);
    assert(d1.extent_int(0) >= n); assert(d2.extent_int(0) >= m);
    assert(x.extent_int(0) >= n && x.extent_int(1) >= nlev);
    assert(y.extent_int(0) >= m && y.extent_int(1) >= nlev);
    using Kokkos::parallel_for;
    const auto ttrn = Kokkos::TeamThreadRange(team, n);
    const auto ttrm = Kokkos::TeamThreadRange(team, m);
    const auto tvr = Kokkos::ThreadVectorRange(team, nlev);
    loop_ik(ttrn, tvr, [&] (int i, int k) { w(i,k) = x(i,k) * d1(i); });
    team.team_barrier();
    parallel_for( ttrm,   [&] (const int i) {
      parallel_for(tvr,   [&] (const int k) { y(i,k) = 0; });
      for (int j = 0; j < n; ++j)
        parallel_for(tvr, [&] (const int k) { y(i,k) += A(i,j) * w(j,k); });
      parallel_for(tvr,   [&] (const int k) { y(i,k) /= s2 * d2(i); }); });
  }

  // Handle (dof,d) vs (d,dof) index ordering.
  template <bool idx_dof_d> static KOKKOS_INLINE_FUNCTION void
  remapd_idx_order (const int dof, const int d, int& i1, int& i2)
  { if (idx_dof_d) { i1 = dof; i2 = d; } else { i1 = d; i2 = dof; } }
  template <bool idx_dof_d> static KOKKOS_INLINE_FUNCTION int
  remapd_idx_dof (const int n, const int idx) { return idx_dof_d ? idx / 2 : idx % n; }
  template <bool idx_dof_d> static KOKKOS_INLINE_FUNCTION int
  remapd_idx_d   (const int n, const int idx) { return idx_dof_d ? idx % 2 : idx / n; }

  /* Compute (1-based indexing)
         xt(i,d,k) = Dinv(i,:,:) x(i,:,k), i = 1:n
         yt(1:m,d,k) = A xt(1:n,d,k)
         y(i,d,k) = D(i,:,:) (yt(i,:,k)/s2),    i = 1:n
     for k = 1:nlev. Sizes are min; a dim can have larger size.
         A m by n,
         Dinv n by 2 by 2, D m by 2 by 2
         x (n,2) by nlev, w (n,2) by nlev with x indexing, y (m,2) by nlev
     Permitted aliases:
         w = x
     If x_idx_dof_d, then x is indexed as x(i,d,k) and y as y(d,i,k); else the
     opposite.
   */
  template <bool x_idx_dof_d,
            typename AT, typename DinvT, typename DT,
            typename XT, typename WT, typename YT>
  static KOKKOS_FUNCTION void
  remapd (const MT& team, const int m, const int n, const int nlev, const AT& A,
          const DinvT& Dinv, const Real s2, const DT& D,
          const XT& x, const WT& w, const YT& y) {
    using Kokkos::parallel_for;
    const auto ttrn  = Kokkos::TeamThreadRange(team,   n);
    const auto ttrm  = Kokkos::TeamThreadRange(team,   m);
    const auto ttr2m = Kokkos::TeamThreadRange(team, 2*m);
    const auto tvr = Kokkos::ThreadVectorRange(team, nlev);
    parallel_for(ttrn, [&] (const int i) {
      // This impl permits w to alias x. The alternative is to use twice as many
      // threads but w can't alias x.
      int i11, i12; remapd_idx_order<x_idx_dof_d>(i, 0, i11, i12);
      int i21, i22; remapd_idx_order<x_idx_dof_d>(i, 1, i21, i22);
      parallel_for(tvr, [&] (const int k) {
        const auto x1 = x(i11,i12,k), x2 = x(i21,i22,k);
        w(i11,i12,k) = Dinv(i,0,0)*x1 + Dinv(i,0,1)*x2;
        w(i21,i22,k) = Dinv(i,1,0)*x1 + Dinv(i,1,1)*x2;
      });
    });
    team.team_barrier();
    parallel_for(ttr2m, [&] (const int idx) {
      const int i = remapd_idx_dof<!x_idx_dof_d>(m, idx);
      const int d = remapd_idx_d  <!x_idx_dof_d>(m, idx);
      int yi1, yi2; remapd_idx_order<!x_idx_dof_d>(i, d, yi1, yi2);
      parallel_for(tvr, [&] (const int k) { y(yi1,yi2,k) = 0; });
      for (int j = 0; j < n; ++j) {
        int xj1, xj2; remapd_idx_order<x_idx_dof_d>(j, d, xj1, xj2);
        parallel_for(tvr, [&] (const int k) { y(yi1,yi2,k) += A(i,j) * w(xj1,xj2,k); });
      }
      parallel_for(tvr, [&] (const int k) { y(yi1,yi2,k) /= s2; });
    });
    team.team_barrier();
    parallel_for(ttrm, [&] (const int i) {
      // This impl avoids a work slot having y's structure; the alternative
      // using twice as many threads requires an extra work slot.
      int i11, i12; remapd_idx_order<!x_idx_dof_d>(i, 0, i11, i12);
      int i21, i22; remapd_idx_order<!x_idx_dof_d>(i, 1, i21, i22);
      parallel_for(tvr, [&] (const int k) {
        const auto y1 = y(i11,i12,k), y2 = y(i21,i22,k);
        y(i11,i12,k) = D(i,0,0)*y1 + D(i,0,1)*y2;
        y(i21,i22,k) = D(i,1,0)*y1 + D(i,1,1)*y2;
      });
    });
  }

  /* CAAS as described in Alg 3.1 of doi:10.1137/18M1165414. q is a mixing
     ratio. Solve
        min_q* norm(dp q - dp q*, 1)
         st    spheremp'(dp q*) = spheremp'(dp q)
               qmin < q* < qmax
     where spheremp = s geo. Operate on DOF indices 0:n-1 and level pack indices
     0:nlev-1.
   */
  template <typename CR1, typename V1, typename CV2, typename V2, typename VQ>
  static KOKKOS_FUNCTION void
  limiter_clip_and_sum (const MT& team, const int n, const int nlev,
                        const Real s, const CR1& geo, const V1& qmin, const V1& qmax,
                        const CV2& dp, const V2& wrk, const VQ& q) {
    assert(geo.extent_int(0) >= n);
    assert(qmin.extent_int(0) >= nlev); assert(qmax.extent_int(0) >= nlev);
    assert(dp .extent_int(0) >= n && dp .extent_int(1) >= nlev);
    assert(wrk.extent_int(0) >= n && wrk.extent_int(1) >= nlev);
    assert(q  .extent_int(0) >= n && q  .extent_int(1) >= nlev);
    static_assert(Scalar::vector_length == packn, "vector_length == packn");
    const auto f = [&] (const int k) {
      auto& c = wrk;
      { // In the case of an infeasible problem, prefer to conserve mass and
        // violate a bound.
        Scalar mass(0), qmass(0);
        for (int i = 0; i < n; ++i) {
          c(i,k) = (s*geo(i))*dp(i,k);
          mass  += c(i,k);
          qmass += c(i,k)*q(i,k);
        }
        VECTOR_SIMD_LOOP for (int s = 0; s < packn; ++s)
          if (qmass[s] < qmin(k)[s]*mass[s])
            qmin(k)[s] = qmass[s]/mass[s];
        VECTOR_SIMD_LOOP for (int s = 0; s < packn; ++s)
          if (qmass[s] > qmax(k)[s]*mass[s])
            qmax(k)[s] = qmass[s]/mass[s];
      }

      Scalar addmass(0);
      bool modified[packn] = {0};
      // Clip.
      for (int i = 0; i < n; ++i)
        VECTOR_SIMD_LOOP for (int s = 0; s < packn; ++s) {
          auto& x = q(i,k)[s];
          const auto xmin = qmin(k)[s];
          const auto xmax = qmax(k)[s];
          if (x > xmax) {
            modified[s] = true;
            addmass[s] += (x - xmax)*c(i,k)[s];
            x = xmax;
          } else if (x < xmin) {
            modified[s] = true;
            addmass[s] += (x - xmin)*c(i,k)[s];
            x = xmin;
          }
        }

      {
        // Compute weights normalization.
        Scalar den(0);
        for (int i = 0; i < n; ++i)
          VECTOR_SIMD_LOOP for (int s = 0; s < packn; ++s)
            if (modified[s]) {
              if (addmass[s] > 0)
                den[s] += (qmax(k)[s] - q(i,k)[s])*c(i,k)[s];
              else
                den[s] += (q(i,k)[s] - qmin(k)[s])*c(i,k)[s];
            }
        // Redistribute mass.
        for (int i = 0; i < n; ++i)
          VECTOR_SIMD_LOOP for (int s = 0; s < packn; ++s)
            if (modified[s] && den[s] > 0) {
              auto& x = q(i,k)[s];
              const auto v = addmass[s] > 0 ? qmax(k)[s] - x : x - qmin(k)[s];
              x += addmass[s]*(v/den[s]);
            }
      }
    };
    team_parallel_for_with_linear_index(team, nlev, f);
  }

  template <typename CR1, typename VW, typename VQ>
  static KOKKOS_FUNCTION void
  limiter_clip_and_sum_real1 (const MT& team, const int n, const Real s, const CR1& geo,
                              Real& qmin, Real& qmax, const VW& wrk, const VQ& q) {
    const auto f = [&] (const int) {
      auto& c = wrk;
      Real mass = 0, qmass = 0;
      for (int i = 0; i < n; ++i) {
        c(i) = s*geo(i);
        mass  += c(i);
        qmass += c(i)*q(i);
      }
      if (qmass < qmin*mass) qmin = qmass/mass;
      if (qmass > qmax*mass) qmax = qmass/mass;

      Real addmass = 0;
      bool modified = false;
      // Clip.
      for (int i = 0; i < n; ++i) {
        auto& x = q(i);
        const auto xmin = qmin;
        const auto xmax = qmax;
        if (x > xmax) {
          modified = true;
          addmass += (x - xmax)*c(i);
          x = xmax;
        } else if (x < xmin) {
          modified = true;
          addmass += (x - xmin)*c(i);
          x = xmin;
        }
      }

      if ( ! modified) return;

      // Compute weights normalization.
      Real den = 0;
      if (addmass > 0) { for (int i = 0; i < n; ++i) den += (qmax - q(i))*c(i); }
      else             { for (int i = 0; i < n; ++i) den += (q(i) - qmin)*c(i); }
      if (den == 0) return;
      // Redistribute mass.
      for (int i = 0; i < n; ++i) {
        const auto v = addmass > 0 ? qmax - q(i) : q(i) - qmin;
        q(i) += addmass*(v/den);
      }
    };
    team_parallel_for_with_linear_index(team, 1, f);
  }

  template <typename View> static KOKKOS_INLINE_FUNCTION
  Real* pack2real (const View& v) { return &(*v.data())[0]; }
  template <typename View> static KOKKOS_INLINE_FUNCTION
  const Real* cpack2real (const View& v) { return &(*v.data())[0]; }
  template <typename View> static KOKKOS_INLINE_FUNCTION
  Scalar* real2pack (const View& v) { return reinterpret_cast<Scalar*>(v.data()); }
  template <typename View> static KOKKOS_INLINE_FUNCTION
  const Scalar* creal2pack (const View& v) { return reinterpret_cast<const Scalar*>(v.data()); }

  // ||4 f(k) on k = 0:niter-1.
  template <typename Fn> KOKKOS_INLINE_FUNCTION static void
  team_parallel_for_with_linear_index (const MT& team, const int niter, const Fn& fn) {
    using Kokkos::parallel_for;
    // 2D -> 1D thread index space. Need to make vector dim the faster for
    // coalesced memory access on the GPU. Kokkos doesn't expose the number of
    // threads in a team, so we have to go to the lower-level API here.
    const int nthr_per_team =
#if defined __CUDA_ARCH__ || defined __HIP_DEVICE_COMPILE__
      blockDim.x,
#else
      1,
#endif
      team_niter = (niter + nthr_per_team - 1)/nthr_per_team;
    assert(OnGpu<ExecSpace>::value || nthr_per_team == 1);
    parallel_for(Kokkos::TeamThreadRange(team, team_niter), [&] (const int team_idx) {
      parallel_for(Kokkos::ThreadVectorRange(team, nthr_per_team), [&] (const int vec_idx) {
        const int k = team_idx*nthr_per_team + vec_idx;
        if (k >= niter) return;
        fn(k);
      }); });
  }

  template <typename TTR, typename TVR, typename Fn>
  static KOKKOS_INLINE_FUNCTION void
  loop_ik (const TTR& ttr, const TVR& tvr, const Fn& fn) {
    using Kokkos::parallel_for;
    parallel_for(ttr, [&] (const int i) { parallel_for(tvr, [&] (const int k) { fn(i,k); }); });
  }
};

} // namespace Homme

#endif // HOMMEXX_GLLFVREMAP_IMPL_HPP
