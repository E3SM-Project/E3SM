/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_DIRK_FUNCTOR_IMPL_HPP
#define HOMMEXX_DIRK_FUNCTOR_IMPL_HPP

#include "Types.hpp"
#include "EquationOfState.hpp"
#include "FunctorsBuffersManager.hpp"
#include "Elements.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "PhysicalConstants.hpp"
#include "ElementOps.hpp"
#include "profiling.hpp"
#include "ErrorDefs.hpp"
#include "utilities/scream_tridiag.hpp"

#include <cassert>

namespace Homme {

struct DirkFunctorImpl {
  enum : int { packn = VECTOR_SIZE };
  enum : int { scaln = NP*NP };
  enum : int { npack = (scaln + packn - 1)/packn };
  enum : int { max_num_lev_pack = NUM_LEV_P };
  enum : int { num_lev_aligned = max_num_lev_pack*packn };
  enum : int { num_phys_lev = NUM_PHYSICAL_LEV };
  enum : int { num_work = 12 };
  enum : bool { calc_initial_guess_in_newton_kernel = false };

  enum : int {
#ifdef HOMMEXX_BFB_TESTING
    default_bfb_solver = true
#else
    default_bfb_solver = false
#endif
  };

  static_assert(num_lev_aligned >= 3,
                "We use wrk(0:2,:) and so need num_lev_aligned >= 3");

  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MT = typename TeamPolicy::member_type;

  using Work
    = Kokkos::View<Scalar*[num_work][num_lev_aligned][npack],
                   Kokkos::LayoutRight, ExecSpace>;
  using WorkSlot
    = Kokkos::View<Scalar           [num_lev_aligned][npack],
                   Kokkos::LayoutRight, ExecSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using ConstWorkSlot
    = Kokkos::View<const Scalar     [num_lev_aligned][npack],
                   Kokkos::LayoutRight, ExecSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using LinearSystem
    = Kokkos::View<Scalar*[4][num_phys_lev][npack],
                   Kokkos::LayoutRight, ExecSpace>;
  using LinearSystemSlot
    = Kokkos::View<Scalar    [num_phys_lev][npack],
                   Kokkos::LayoutRight, ExecSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  KOKKOS_INLINE_FUNCTION
  static WorkSlot get_work_slot (const Work& w, const int& wi, const int& si) {
    using Kokkos::subview;
    using Kokkos::ALL;
    const auto a = ALL();
    return subview(w, wi, si, a, a);
  }

  KOKKOS_INLINE_FUNCTION
  static LinearSystemSlot get_ls_slot (const LinearSystem& w, const int& wi,
                                       const int& si) {
    using Kokkos::subview;
    using Kokkos::ALL;
    const auto a = ALL();
    return subview(w, wi, si, a, a);
  }

  Work m_work;
  LinearSystem m_ls;
  TeamPolicy m_policy, m_ig_policy;
  TeamUtils<ExecSpace> m_tu, m_tu_ig;
  int nteam;

  KOKKOS_INLINE_FUNCTION
  size_t shmem_size (const int team_size) const {
    return KernelVariables::shmem_size(team_size);
  }

  DirkFunctorImpl (const int nelem)
    : m_policy(1,1,1), m_ig_policy(1,1,1), m_tu(m_policy), m_tu_ig(m_ig_policy) // throwaway settings
  {
    init(nelem);
  }

  void init (const int nelem) {
    if (OnGpu<ExecSpace>::value) {
      ThreadPreferences tp;
      tp.max_threads_usable = NUM_PHYSICAL_LEV;
      tp.max_vectors_usable = NP*NP;
      tp.prefer_threads = false;
      const auto p = DefaultThreadsDistribution<ExecSpace>
        ::team_num_threads_vectors(nelem, tp);
      const auto
        nhwthr = p.first*p.second,
        nvec = std::min(NP*NP, nhwthr),
        nthr = nhwthr/nvec;
      m_policy = TeamPolicy(nelem, nthr, nvec);
    } else {
      ThreadPreferences tp;
      tp.max_threads_usable = NUM_PHYSICAL_LEV;
      tp.max_vectors_usable = 1;
      tp.prefer_threads = true;
      const auto p = DefaultThreadsDistribution<ExecSpace>
        ::team_num_threads_vectors(nelem, tp);
      m_policy = TeamPolicy(nelem, p.first, 1);
    }
    m_tu = TeamUtils<ExecSpace>(m_policy);
    nteam = std::min(nelem, m_tu.get_num_concurrent_teams());
    m_ig_policy = Homme::get_default_team_policy<ExecSpace>(nelem);
    m_tu_ig = TeamUtils<ExecSpace>(m_ig_policy);
  }

  int requested_buffer_size () const {
    // FunctorsBuffersManager wants the size in terms of sizeof(Real).
    return (Work::shmem_size(nteam) + LinearSystem::shmem_size(nteam))/sizeof(Real);
  }

  void init_buffers (const FunctorsBuffersManager& fbm) {
    Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());
    m_work = Work(mem, nteam);
    mem += Work::shmem_size(nteam)/sizeof(Scalar);
    m_ls = LinearSystem(mem, nteam);
  }

  void run (int nm1, Real alphadt_nm1, int n0, Real alphadt_n0, int np1, Real dt2,
            const Elements& e, const HybridVCoord& hvcoord,
            const bool bfb_solver = default_bfb_solver) {
    if ( ! calc_initial_guess_in_newton_kernel) {
      run_initial_guess(np1, e, hvcoord);
      Kokkos::fence();
    }

    run_newton(nm1, alphadt_nm1, n0, alphadt_n0, np1, dt2, e, hvcoord, bfb_solver);
    Kokkos::fence();
  }

  // Optimal impl of phi_from_eos for the initial guess. See comments for the
  // function phi_from_eos, below, for discussion. This kernel uses standard
  // Hommexx layout and parallelization approaches to compute the scans
  // optimally.
  void run_initial_guess (int np1, const Elements& e, const HybridVCoord& hvcoord) {
    using Kokkos::subview;
    using Kokkos::parallel_for;

    const auto e_phis = e.m_geometry.m_phis;
    const auto e_vtheta_dp = e.m_state.m_vtheta_dp;
    const auto e_dp3d = e.m_state.m_dp3d;
    const auto e_phinh_i = e.m_derived.m_divdp_proj;

    ElementOps elem_ops;
    elem_ops.init(hvcoord);

    const auto work = m_work;

    const auto toplevel = KOKKOS_LAMBDA (const MT& team) {
      KernelVariables kv(team, m_tu_ig);
      const auto ie = kv.ie;

      const ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV_P]>
        p_i  (reinterpret_cast<Scalar*>(get_work_slot(work, kv.team_idx, 0).data())),
        phi_i(reinterpret_cast<Scalar*>(get_work_slot(work, kv.team_idx, 2).data()));
      const ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV  ]>
        p    (reinterpret_cast<Scalar*>(get_work_slot(work, kv.team_idx, 1).data()));

      const auto f = [&] (const int idx) {
        const int igp = idx / NP, jgp = idx % NP;

        const auto p_i_c = Homme::subview(p_i,igp,jgp);
        const auto p_c = Homme::subview(p,igp,jgp);
        const auto phi_i_c = Homme::subview(phi_i,igp,jgp);

        elem_ops.compute_hydrostatic_p(kv, Homme::subview(e_dp3d,ie,np1,igp,jgp), p_i_c, p_c);
        EquationOfState::compute_phi_i(kv, e_phis(ie,igp,jgp),
                                       Homme::subview(e_vtheta_dp,ie,np1,igp,jgp),
                                       p_c, phi_i_c);

        // Copy phi_i for use in run_newton. Don't need nlevp, since that's
        // always phis.
        const auto g = [&] (const int k) { e_phinh_i(ie,igp,jgp,k) = phi_i_c(k); };
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), g);
      };
      parallel_for(Kokkos::TeamThreadRange(kv.team, NP*NP), f);
    };
    Kokkos::parallel_for(m_ig_policy, toplevel);
  }

  void run_newton (int nm1, Real alphadt_nm1, int n0, Real alphadt_n0, int np1, Real dt2,
                   const Elements& e, const HybridVCoord& hvcoord, const bool bfb_solver) {
    using Kokkos::subview;
    using Kokkos::parallel_for;
    const auto a = Kokkos::ALL();

    const auto grav = PhysicalConstants::g;
    const int nvec = npack;
    const int maxiter = 20;
#ifdef HOMMEXX_BFB_TESTING
    const Real deltatol = 1e-6; // In bfb testing, use coarse tolerance, due to zeroulp calls
#else
    const Real deltatol = 1e-11; // exit if newton increment < deltatol
#endif

    const auto work = m_work;
    const auto ls = m_ls;
    const auto e_w_i = e.m_state.m_w_i;
    const auto e_vtheta_dp = e.m_state.m_vtheta_dp;
    const auto e_phinh_i = e.m_state.m_phinh_i;
    const auto e_dp3d = e.m_state.m_dp3d;
    const auto e_v = e.m_state.m_v;
    const auto e_phis = e.m_geometry.m_phis;
    const auto e_gradphis = e.m_geometry.m_gradphis;
    const auto e_initial_guess = e.m_derived.m_divdp_proj;
    const auto hybi = hvcoord.hybrid_bi;

    const auto toplevel = KOKKOS_LAMBDA (const MT& team) {
      KernelVariables kv(team, m_tu);
      const auto ie = kv.ie;
      const int nlev = num_phys_lev;

      const auto
      phi_n0    = get_work_slot(work, kv.team_idx,  0),
      phi_np1   = get_work_slot(work, kv.team_idx,  1),
      dphi      = get_work_slot(work, kv.team_idx,  2),
      w_n0      = get_work_slot(work, kv.team_idx,  3),
      w_np1     = get_work_slot(work, kv.team_idx,  4),
      dpnh_dp_i = get_work_slot(work, kv.team_idx,  5),
      gwh_i     = get_work_slot(work, kv.team_idx,  6),
      dphi_n0   = get_work_slot(work, kv.team_idx,  6), // reuse gwh_i
      vtheta_dp = get_work_slot(work, kv.team_idx,  7),
      dp3d      = get_work_slot(work, kv.team_idx,  8),
      pnh       = get_work_slot(work, kv.team_idx,  9),
      wrk       = get_work_slot(work, kv.team_idx, 10),
      xfull     = get_work_slot(work, kv.team_idx, 11);
      const auto
      dl = get_ls_slot(ls, kv.team_idx, 0),
      d  = get_ls_slot(ls, kv.team_idx, 1),
      du = get_ls_slot(ls, kv.team_idx, 2);

      // View of xfull for use in the solver. We want xfull so that we
      // can use the nlevp-1 entry, which we make sure is 0, when convenient.
      LinearSystemSlot x = subview(xfull, Kokkos::pair<int,int>(0,nlev), a);
      xfull(nlev,0)[0] = 0.0;

      const auto transpose4 = [&] (const int nt, const bool transpose_phi_np1 = true) {
        transpose(kv, nlev+1, subview(e_w_i      ,ie,nt,a,a,a), w_np1    );
        transpose(kv, nlev,   subview(e_vtheta_dp,ie,nt,a,a,a), vtheta_dp);
        transpose(kv, nlev,   subview(e_dp3d     ,ie,nt,a,a,a), dp3d     );
        if ( ! transpose_phi_np1) return;
        transpose(kv, nlev+1, subview(e_phinh_i  ,ie,nt,a,a,a), phi_np1  );
      };

      const auto accum_n0 = [&] (const Real dt3, const int nt) {
        // Computing these transposes is inefficient, but doing so lets us use
        // the same pnh_and_exner_from_eos as we use elsewhere. This function is
        // called only in ~1 out of 5 DIRK stage calls, and only outside of the
        // Newton iteration. (Also, tbc, transpose and pnh_and_exner_from_eos
        // are parallel efficient.)
        transpose4(nt);
        calc_gwphis(kv, subview(e_dp3d,ie,nt,a,a,a), subview(e_v,ie,nt,a,a,a,a),
                    subview(e_gradphis,ie,a,a,a), hybi, gwh_i);
        kv.team_barrier();
        loop_ki(kv, nlev, nvec, [&] (int k, int i) {
          dphi(k,i) = phi_np1(k+1,i) - phi_np1(k,i);
        });
        kv.team_barrier();
        pnh_and_exner_from_eos(kv, hvcoord, vtheta_dp, dp3d, dphi, pnh, wrk, dpnh_dp_i);
        kv.team_barrier();
        loop_ki(kv, nlev, nvec, [&] (const int k, const int i) {
          w_n0(k,i) += dt3*grav*(dpnh_dp_i(k,i) - 1);
        });
        loop_ki(kv, nlev, nvec, [&] (int k, int i) {
          phi_n0(k,i) = phi_n0(k,i) + dt3*grav*w_np1(k,i) - dt3*gwh_i(k,i);
        });
      };

      // Compute w_n0, phi_n0.
      transpose(kv, nlev+1, subview(e_phinh_i,ie,np1,a,a,a), phi_n0);
      transpose(kv, nlev+1, subview(e_w_i    ,ie,np1,a,a,a), w_n0  );
      kv.team_barrier();
      // wmax is computed before optional updates to w_n0.
      const auto wmax = calc_wmax(kv, nlev+1, nvec, w_n0);
      // Computed only in some cases.
      if (alphadt_n0 != 0) {
        accum_n0(alphadt_n0, n0);
        kv.team_barrier();
      }
      if (alphadt_nm1 != 0) {
        assert(nm1 >= 0);
        accum_n0(alphadt_nm1, nm1);
        kv.team_barrier();
      }
      // Always computed.
      transpose4(np1, false);
      calc_gwphis(kv, subview(e_dp3d,ie,np1,a,a,a), subview(e_v,ie,np1,a,a,a,a),
                  subview(e_gradphis,ie,a,a,a), hybi, gwh_i);
      kv.team_barrier();
      loop_ki(kv, nlev, nvec, [&] (int k, int i) { phi_n0(k,i) -= dt2*gwh_i(k,i); });

      // Initial guess for phi_np1.
      if (calc_initial_guess_in_newton_kernel) {
        // Use hydrostatic phi.
        phi_from_eos(kv, nlev, nvec, hvcoord, subview(e_phis,ie,a,a), vtheta_dp, dp3d, phi_np1);
      } else {
        // Copy initial guess from where run_initial_guess stashed it.
        transpose(kv, nlev, subview(e_initial_guess,ie,a,a,a), phi_np1);
        loop_ki(kv, 1, nvec, [&] (int, int i) { set_phis(i, subview(e_phis,ie,a,a), phi_np1); });
      }
      kv.team_barrier();
      loop_ki(kv, nlev, nvec, [&] (int k, int i) { dphi(k,i) = phi_np1(k+1,i) - phi_np1(k,i); });
      kv.team_barrier();
      // If any dphi > -g in a column, set it to -g and integrate to get a
      // new initial phi_np1 and w_np1.
      calc_whether_gt_and_set(kv, nlev, nvec, -grav, dphi, wrk);
      kv.team_barrier();
      if (wrk(1,0)[0] == 1) {
        scan_dphi(kv, nlev, nvec, wrk, dphi, phi_np1);
        kv.team_barrier();
      }

      // Initial guess for w_np1.
      loop_ki(kv, nlev, nvec, [&] (int k, int i) { w_np1(k,i) = (phi_np1(k,i) - phi_n0(k,i))/(dt2*grav); });

      loop_ki(kv, nlev, nvec, [&] (int k, int i) { dphi_n0(k,i) = phi_n0(k+1,i) - phi_n0(k,i); });

      int it = 0;
      Real deltaerr;
      for (; it < maxiter; ++it) { // Newton iteration
        pnh_and_exner_from_eos(kv, hvcoord, vtheta_dp, dp3d, dphi, pnh, wrk, dpnh_dp_i);
        kv.team_barrier();
        loop_ki(kv, nlev, nvec, [&] (const int k, const int i) {
          x(k,i) = -(w_np1(k,i) - (w_n0(k,i) + grav*dt2*(dpnh_dp_i(k,i) - 1))); // -residual
        });

        calc_jacobian(kv, dt2, dp3d, dphi, pnh, dl, d, du);
        kv.team_barrier();
        if (bfb_solver) solvebfb(kv, dl, d, du, x); else solve(kv, dl, d, du, x);
        kv.team_barrier();

        loop_ki(kv, 1, nvec, [&] (int k, int i) { wrk(2,i) = 1; });
        kv.team_barrier();
        for (int nsafe = 0; nsafe < 2; ++nsafe) {
          loop_ki(kv, nlev-1, nvec, [&] (int k, int i) {
            dphi(k,i) = dphi_n0(k,i) + dt2*grav*(         (w_np1(k+1,i) - w_np1(k,i)) +
                                                 wrk(2,i)*(    x(k+1,i) -     x(k,i)));
          });
          loop_ki(kv, 1, nvec, [&] (int, int i) {
            const auto k = nlev-1;
            dphi(k,i) = dphi_n0(k,i) - dt2*grav*(w_np1(k,i) + wrk(2,i)*x(k,i));
          });
          kv.team_barrier();
          calc_whether_ge(kv, nlev, nvec, 0, dphi, wrk);
          kv.team_barrier();
          if (wrk(1,0)[0] == 0) break;
          calc_step_size(kv, nlev, nvec, grav, dt2, dphi_n0, w_np1, x, wrk);
          kv.team_barrier();
        }
        kv.team_barrier();

        loop_ki(kv, nlev, nvec, [&] (int k, int i) { w_np1(k,i) += wrk(2,i)*x(k,i); });

        if (exit_on_step(kv, nlev, nvec, wmax, deltatol, x, deltaerr)) break;
      } // Newton iteration
      kv.team_barrier();

      if (it>=maxiter) {
        printf ("[DIRK] WARNING! Newton reached max iteration count, with deltaerr = %3.17f\n",deltaerr);
      }

      // Update phi_np1.
      loop_ki(kv, nlev, nvec, [&] (int k, int i) { phi_np1(k,i) = phi_n0(k,i) + dt2*grav*w_np1(k,i); });

      kv.team_barrier();
      transpose(kv, nlev+1, phi_np1, subview(e_phinh_i,ie,np1,a,a,a));
      transpose(kv, nlev+1, w_np1,   subview(e_w_i    ,ie,np1,a,a,a));
    };

    Kokkos::parallel_for(m_policy, toplevel);
  }

  template <typename Fn>
  KOKKOS_INLINE_FUNCTION
  static void loop_ki (const KernelVariables& kv, const int klim,
                       const int ilim, const Fn& g) {
    using Kokkos::parallel_for;
    using Kokkos::TeamThreadRange;
    using Kokkos::ThreadVectorRange;

    if (OnGpu<ExecSpace>::value) {
      const auto tr = TeamThreadRange  (kv.team, klim);
      const auto vr = ThreadVectorRange(kv.team, ilim);
      const auto f = [&] (const int k) {
        const auto h = [&] (const int i) { g(k,i); };
        parallel_for(vr, h);
      };
      parallel_for(tr, f);
    } else if (kv.team.team_size() == 1) {
      for (int k = 0; k < klim; ++k)
        for (int i = 0; i < ilim; ++i)
          g(k,i);
    } else {
      const auto tr = TeamThreadRange  (kv.team, klim);
      const auto f = [&] (const int k) {
        for (int i = 0; i < ilim; ++i)
          g(k,i);
      };
      parallel_for(tr, f);
    }
  }

  // Format of rest of Hxx -> DIRK Newton iteration format.
  template <typename View>
  KOKKOS_INLINE_FUNCTION
  static void transpose (const KernelVariables& kv, const int nlev,
                         const View& src, const WorkSlot& dst,
                         typename std::enable_if<View::rank == 3>::type* = 0) {
    assert(src.extent_int(2)*packn >= nlev);
    assert(src.extent_int(0) == NP && src.extent_int(1) == NP);
    const auto f = [&] (const int k) {
      const auto
      pk = k / packn,
      sk = k % packn;
      const auto g = [&] (const int i) {
        const auto gk0 = packn*i;
        for (int s = 0; s < packn; ++s) {
          const auto
          gk = gk0 + s,
          gi = gk / NP,
          gj = gk % NP;
          if (scaln % packn != 0 && // try to compile out this conditional when possible
              gk >= scaln) break;
          dst(k,i)[s] = src(gi,gj,pk)[sk];
        }
      };
      const int n = npack;
      const auto p = Kokkos::ThreadVectorRange(kv.team, n);
      Kokkos::parallel_for(p, g);
    };    
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, nlev), f);
  }

  // DIRK Newton iteration format -> format of rest of Hxx.
  template <typename View>
  KOKKOS_INLINE_FUNCTION
  static void transpose (const KernelVariables& kv, const int nlev,
                         const WorkSlot& src, const View& dst,
                         typename std::enable_if<View::rank == 3>::type* = 0) {
    assert(dst.extent_int(2)*packn >= nlev);
    assert(dst.extent_int(0) == NP && dst.extent_int(1) == NP);
    const auto f = [&] (const int idx) {
      const auto
      gi = idx / NP,
      gj = idx % NP,
      pi = idx / packn,
      si = idx % packn;
      const auto g = [&] (const int pk) {
        const auto k0 = pk*packn;
        // If there is a remainder at the end, we nonetheless transfer these
        // unused data to avoid a conditional and runtime loop limit.
        for (int sk = 0; sk < packn; ++sk)
          dst(gi,gj,pk)[sk] = src(k0+sk,pi)[si];
      };
      const auto p = Kokkos::ThreadVectorRange(kv.team, dst.extent_int(2));
      Kokkos::parallel_for(p, g);
    };    
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP*NP), f);
  }

  // Compute a vertical velocity induced by surface topography:
  //     wh_i = ubar grad phis
  template <typename R, typename Rv, typename Rgphis, typename Rhybi, typename W>
  KOKKOS_INLINE_FUNCTION
  static void calc_gwphis (
    const KernelVariables& kv,
    // All in arrays are in Hxx format.
    const R& dp3d, const Rv& v, const Rgphis& gradphis, const Rhybi& hybi,
    // Out array is in DIRK format. gwh_i(nlevp,:) is not written since it is
    // not used.
    const W& gwh_i,
    const int nlev = NUM_PHYSICAL_LEV)
  {
    using Kokkos::parallel_for;
    using Kokkos::TeamThreadRange;
    using Kokkos::ThreadVectorRange;

    const int n = npack;
    const auto pv = ThreadVectorRange(kv.team, n);

    const auto f1 = [&] (const int) {
      const auto k0 = [&] (const int i) { gwh_i(0,i) = 0; };
      parallel_for(pv, k0);
    };
    parallel_for(TeamThreadRange(kv.team, 1), f1);

    const auto f2 = [&] (const int km1) {
      const auto
      k = km1 + 1,
      pk = k / packn,
      sk = k % packn,
      pkm1 = (k-1) / packn,
      skm1 = (k-1) % packn;
      const auto g = [&] (const int i) {
        Scalar dp3dk, dp3dkm1, v1k, v2k, v1km1, v2km1, gphis1, gphis2;
        for (int s = 0; s < packn; ++s) {
          const auto
          idx = packn*i + s,
          gi = idx / NP,
          gj = idx % NP;
          if (scaln % packn != 0 && idx >= scaln) break;
          dp3dkm1[s] = dp3d(gi,gj,pkm1)[skm1];
          dp3dk  [s] = dp3d(gi,gj,pk  )[sk];
          v1km1  [s] = v (0,gi,gj,pkm1)[skm1];
          v2km1  [s] = v (1,gi,gj,pkm1)[skm1];
          v1k    [s] = v (0,gi,gj,pk  )[sk];
          v2k    [s] = v (1,gi,gj,pk  )[sk];
          gphis1 [s] = gradphis(0,gi,gj);
          gphis2 [s] = gradphis(1,gi,gj);
        }
        const auto den = dp3dkm1 + dp3dk;
        const auto v1_i = (dp3dk*v1k + dp3dkm1*v1km1) / den;
        const auto v2_i = (dp3dk*v2k + dp3dkm1*v2km1) / den;
        gwh_i(k,i) = (v1_i*gphis1 + v2_i*gphis2)*hybi(k);
      };
      parallel_for(pv, g);
    };
    parallel_for(TeamThreadRange(kv.team, nlev-1), f2);
  }

  template <typename R, typename W, typename Wi>
  KOKKOS_INLINE_FUNCTION
  static void pnh_and_exner_from_eos (
    const KernelVariables& kv, const HybridVCoord& hvcoord,
    // All arrays are in DIRK format.
    const R& vtheta_dp, const R& dp3d, const R& dphi,
    // exner is workspace. dpnh_dp_i(nlevp,:) is not computed.
    const W& pnh, const W& exner, const Wi& dpnh_dp_i,
    const int nlev = NUM_PHYSICAL_LEV)
  {
    using Kokkos::parallel_for;

    const int n = npack;
    const auto pv = Kokkos::ThreadVectorRange(kv.team, n);

    // Compute pnh(1:nlev,:). pnh(nlevp,:) is not needed.
    const auto f1 = [&] (const int k) {
      const auto g = [&] (const int i) {
        EquationOfState::compute_pnh_and_exner(
          vtheta_dp(k,i), dphi(k,i), pnh(k,i), exner(k,i));
      };
      parallel_for(pv, g);
    };
    parallel_for(Kokkos::TeamThreadRange(kv.team, nlev), f1);

    // Compute dpnh_dp_i(1:nlev,:), d(pnh)/d(pi) at interfaces. Use one-sided
    // differences at boundaries. Do not compute dpnh_dp_i(nlevp,:).
    kv.team_barrier(); // wait for pnh
    const auto f2 = [&] (const int) {
      const auto k0 = [&] (const int i) {
        const auto pnh_i_0 = hvcoord.hybrid_ai0*hvcoord.ps0; // hydrostatic ptop
        dpnh_dp_i(0,i) = 2*(pnh(0,i) - pnh_i_0)/dp3d(0,i);
      };
      parallel_for(pv, k0);
    };
    parallel_for(Kokkos::TeamThreadRange(kv.team, 1), f2);
    const auto f3 = [&] (const int km1) {
      const auto k = km1 + 1;
      const auto kr = [&] (const int i) {
        dpnh_dp_i(k,i) = ((pnh(k,i) - pnh(k-1,i))/
                          ((dp3d(k-1,i) + dp3d(k,i))/2));
      };
      parallel_for(pv, kr);
    };
    parallel_for(Kokkos::TeamThreadRange(kv.team, nlev-1), f3);
  }

  // Suboptimal impl of the initial guess that uses the policy native to the
  // DIRK Newton iteration. It is suboptimal because there are two scans, and
  // both the memory layout and the parallelization are suboptimal for these.
  // run_initial_guess, above, runs a separate kernel with optimal layout and
  // parallelization.
  template <typename Rphis, typename R, typename W>
  KOKKOS_INLINE_FUNCTION static void
  phi_from_eos (const KernelVariables& kv, const int nlev, const int nvec,
                const HybridVCoord& hvcoord, const Rphis& phis, const R& vtheta_dp, const R& dp,
                // phi_i on output
                const W& wrk)
  {
    // Scan to compute pressure.
    loop_ki(kv, 1, nvec, [&] (int, int i) {
      wrk(0,i) = hvcoord.hybrid_ai0*hvcoord.ps0;
      for (int k = 0; k < nlev; ++k) {
        wrk(k+1,i) = wrk(k,i) + dp(k,i);      // p_i at interface k+1
        wrk(k,i) = (wrk(k+1,i) + wrk(k,i))/2; // p_m at midpoint  k
      }
    });
    kv.team_barrier();
    // Do most of the flops.
    loop_ki(kv, nlev, nvec, [&] (int k, int i) {
      wrk(k,i) = EquationOfState::compute_dphi(vtheta_dp(k,i), wrk(k,i)); // dphi
    });
    kv.team_barrier();
    // Scan to compute phi_i.
    loop_ki(kv, 1, nvec, [&] (int, int i) {
      set_phis(i, phis, wrk);
      for (int k = nlev-1; k >= 0; --k)
        wrk(k,i) = wrk(k+1,i) + wrk(k,i); // phi_i below + dphi
    });
  }

  template <typename Rphis, typename W>
  KOKKOS_INLINE_FUNCTION static void
  set_phis (const int i, const Rphis& phis, const W& phi_i) {
    for (int s = 0; s < packn; ++s) {
      const int idx = i*packn + s, gi = idx / NP, gj = idx % NP;
      if (scaln % packn != 0 && idx >= scaln) break;
      phi_i(num_phys_lev,i)[s] = phis(gi,gj);
    }    
  }

  KOKKOS_INLINE_FUNCTION
  static Real calc_wmax (const KernelVariables& kv, const int nlev, const int nvec,
                         const WorkSlot& w) {
    using Kokkos::parallel_reduce;
    using Kokkos::TeamThreadRange;
    using Kokkos::ThreadVectorRange;

    const auto f = [&] (int k, Real& maxval) {
      const auto g = [&] (int i, Real& lmaxval) {
        const auto v = w(k,i);
        for (int s = 0; s < packn; ++s) {
          if (scaln % packn != 0 && i*packn + s >= scaln) break;
          lmaxval = max(lmaxval, std::abs(v[s]));
        }
      };
      Real lmaxval;
      const auto vr = ThreadVectorRange(kv.team, nvec);
      parallel_reduce(vr, g, Kokkos::Max<Real>(lmaxval));
      maxval = max(maxval, lmaxval); // benign write race
    };
    Real wmax;
    const auto tr = TeamThreadRange(kv.team, nlev);
    parallel_reduce(tr, f, Kokkos::Max<Real>(wmax));
    return max(1.0, wmax);
  }

  KOKKOS_INLINE_FUNCTION
  static bool exit_on_step (const KernelVariables& kv, const int nlev, const int nvec,
                            const Real& wmax, const Real& deltatol,
                            const LinearSystemSlot& x, Real& deltaerr) {
    using Kokkos::parallel_reduce;
    using Kokkos::TeamThreadRange;
    using Kokkos::ThreadVectorRange;

    // deltaerr = maxval(abs(x) / wmax
    const auto f = [&] (int k, Real& maxval) {
      const auto g = [&] (int i, Real& lmaxval) {
        const auto v = x(k,i);
        for (int s = 0; s < packn; ++s) {
          if (scaln % packn != 0 && i*packn + s >= scaln) break;
          lmaxval = max(lmaxval, std::abs(v[s]));
        }
      };
      Real lmaxval;
      const auto vr = ThreadVectorRange(kv.team, nvec);
      parallel_reduce(vr, g, Kokkos::Max<Real>(lmaxval));
      maxval = max(maxval, lmaxval); // benign write race
    };
    const auto tr = TeamThreadRange(kv.team, nlev);
    parallel_reduce(tr, f, Kokkos::Max<Real>(deltaerr));
    return deltaerr/wmax < deltatol;
  }

  /* Compute Jacobian of F(phi) = sum(dphi) + const + (dt*g)^2 *(1-dp/dpi)
     column wise with respect to phi. Form the tridiagonal analytical Jacobian J
     to solve J * x = -f.

     This code will need to change when the equation of state is changed.
  */
  template <typename R, typename W>
  KOKKOS_INLINE_FUNCTION
  static void calc_jacobian (const KernelVariables& kv, const Real& dt2,
                             // All arrays are in DIRK format.
                             const R& dp3d, const R& dphi, const R& pnh,
                             const W& dl, const W& d, const W& du,
                             const int nlev = NUM_PHYSICAL_LEV) {
    using Kokkos::parallel_for;
    using Scalar = typename R::value_type;

    const int n = npack;
    const auto pv = Kokkos::ThreadVectorRange(kv.team, n);
    const auto pt1 = Kokkos::TeamThreadRange(kv.team, 1);

    const Real a = square(dt2*PhysicalConstants::g)/(1 - PhysicalConstants::kappa);

    const auto f1 = [&] (const int) {
      const auto ks = [&] (const int i) { // first Jacobian row
        const int k = 0;
        const auto b = a/dp3d(k,i);
        du(k,i) = 2*b*(pnh(k,i)/dphi(k,i));
        d (k,i) = 1 - du(k,i);
      };
      parallel_for(pv, ks);
    };
    parallel_for(pt1, f1);
    const auto f2 = [&] (const int km1) {
      const auto k = km1 + 1;
      const auto kmid = [&] (const int i) { // middle Jacobian rows
        const auto b = 2*a/(dp3d(k-1,i) + dp3d(k,i));
        dl(k,i) = b*(pnh(k-1,i)/dphi(k-1,i));
        du(k,i) = b*(pnh(k  ,i)/dphi(k  ,i));
        // In all rows k,
        //     dl <= 0, du <= 0,
        // and thus
        //     d = 1 + |dl| + |du| > |dl| + |du|,
        // making this Jacobian matrix strictly diagonally dominant. Thus, we
        // need not pivot when factorizing the matrix.
        d (k,i) = 1 - dl(k,i) - du(k,i);
      };
      parallel_for(pv, kmid);
    };
    parallel_for(Kokkos::TeamThreadRange(kv.team, nlev-2), f2);
    const auto f3 = [&] (const int) {
      const auto ke = [&] (const int i) { // last Jacobian row
        const int k = nlev-1;
        const auto b = 2*a/(dp3d(k-1,i) + dp3d(k,i));
        dl(k,i) = b*(pnh(k-1,i)/dphi(k-1,i));
        d (k,i) = 1 - dl(k,i) - b*(pnh(k,i)/dphi(k,i));        
      };
      parallel_for(pv, ke);
    };
    parallel_for(pt1, f3);
  }

  template <typename W>
  KOKKOS_INLINE_FUNCTION
  static void solve (const KernelVariables& kv,
                     const W& dl, const W& d, const W& du, const W& x) {
    assert(d.extent_int(0) == num_phys_lev);
    if (OnGpu<ExecSpace>::value)
      scream::tridiag::cr(kv.team, dl, d, du, x);
    else {
      const auto f = [&] () { scream::tridiag::thomas(dl, d, du, x); };
      Kokkos::single(Kokkos::PerTeam(kv.team), f);
    }
  }

  template <typename W>
  KOKKOS_INLINE_FUNCTION
  static void solvebfb (const KernelVariables& kv,
                        const W& dl, const W& d, const W& du, const W& x) {
    assert(d.extent_int(0) == num_phys_lev);
    scream::tridiag::bfb(kv.team, dl, d, du, x);
  }

  // Determine a step length 0 < alpha <= 1.
  KOKKOS_INLINE_FUNCTION static void
  calc_step_size (const KernelVariables& kv, const int nlev, const int nvec,
                  const Real& grav, const Real& dt2,
                  const WorkSlot& dphi_n0, const WorkSlot& w_np1, const LinearSystemSlot& x,
                  // On input, wrk(0,i)[s] is 1 if the step length should be
                  // smaller. On output, alpha is in wrk(2,:).
                  const WorkSlot& wrk) {
    using Kokkos::parallel_reduce;
    using Kokkos::parallel_for;
    using Kokkos::TeamThreadRange;
    using Kokkos::ThreadVectorRange;
    // Set all alpha_k to 1 except row 0. Include row nlev for use as a flag.
    loop_ki(kv, nlev, nvec, [&] (int k, int i) { wrk(k+1,i) = 1; });
    kv.team_barrier();
    loop_ki(kv, nlev, nvec, [&] (int k, int i) {
      for (int s = 0; s < packn; ++s) {
        if (scaln % packn != 0 && i*packn + s >= scaln) break;
        if (wrk(0,i)[s] == 0) {
          // Indicate this column is already good.
          wrk(nlev,i)[s] = 0;
          continue;
        }
        Real dx, dw;
        if (k < nlev-1) {
          dx =     x(k+1,i)[s] -     x(k,i)[s];
          dw = w_np1(k+1,i)[s] - w_np1(k,i)[s];
        } else {
          dx = -    x(k,i)[s];
          dw = -w_np1(k,i)[s];          
        }
        if (dx != 0) {
          // Step length at which dphi(k,i)[s] would = 0.
          const Real alpha = -(dphi_n0(k,i)[s] + dt2*grav*dw)/(dt2*grav*dx);
          // A negative step is irrelevant.
          if (alpha >= 0) wrk(k,i)[s] = alpha;
        }
      }
    });
    kv.team_barrier();
    // Find minimum alpha in each column. This calculation is performed rarely,
    // so we can use the following suboptimal reduction arising from the
    // Newton-loop team policy without essentially any performance loss. The
    // code itself is fine, but the current policy has too many threads and too
    // few vector lanes on GPU compared with the optimal policy.
    const auto f = [&] (int idx) {
      const int i = idx / packn, s = idx % packn;
      if (wrk(nlev,i)[s] == 0) return;
      const auto g = [&] (int k, Real& lalpha) { lalpha = min(lalpha, wrk(k,i)[s]); };
      Real alpha;
      const auto vr = ThreadVectorRange(kv.team, nlev);
      parallel_reduce(vr, g, Kokkos::Min<Real>(alpha));
      // Step halfway to the distance at which at least one dphi is 0.
      wrk(2,i)[s] = min(1.0, alpha)/2;
    };
    const auto tr = TeamThreadRange(kv.team, static_cast<int>(scaln));
    parallel_for(tr, f);
  }

  // Determine whether any dphi > threshold. Set the entry to threshold.
  KOKKOS_INLINE_FUNCTION static void
  calc_whether_gt_and_set (const KernelVariables& kv, const int nlev, const int nvec,
                           const Real threshold,
                           // dphi > threshold?
                           const WorkSlot& dphi,
                           // On output, wrk(0,:) contains 0 for no, 1 for yes;
                           // if any is yes, then wrk(1,0) is 1, else 0.
                           const WorkSlot& wrk) {
    loop_ki(kv, 2, nvec, [&] (int k, int i) { wrk(k,i) = 0; });
    kv.team_barrier();
    loop_ki(kv, nlev, nvec, [&] (int k, int i) {
      for (int s = 0; s < packn; ++s) {
        if (scaln % packn != 0 && i*packn + s >= scaln) break;
        if (dphi(k,i)[s] > threshold) {
          dphi(k,i)[s] = threshold;
          wrk(0,i)[s] = 1; // benign write race
          wrk(1,0)[0] = 1; // benign write race
        }
      }
    });
  }

  // Determine whether any dphi >= threshold.
  KOKKOS_INLINE_FUNCTION static void
  calc_whether_ge (const KernelVariables& kv, const int nlev, const int nvec,
                   const Real threshold,
                   // dphi >= threshold?
                   const WorkSlot& dphi,
                   // On output, wrk(0,:) contains 0 for no, 1 for yes;
                   // if any is yes, then wrk(1,0) is 1, else 0.
                   const WorkSlot& wrk) {
    loop_ki(kv, 2, nvec, [&] (int k, int i) { wrk(k,i) = 0; });
    kv.team_barrier();
    loop_ki(kv, nlev, nvec, [&] (int k, int i) {
      for (int s = 0; s < packn; ++s) {
        if (scaln % packn != 0 && i*packn + s >= scaln) break;
        if (dphi(k,i)[s] >= threshold) {
          wrk(0,i)[s] = 1; // benign write race
          wrk(1,0)[0] = 1; // benign write race
        }
      }
    });
  }

  KOKKOS_INLINE_FUNCTION static bool
  scan_dphi (const KernelVariables& kv, const int nlev, const int nvec,
             // On input, wrk(0,:) contains 0 for no, 1 for yes
             const WorkSlot& wrk, const WorkSlot& dphi, const WorkSlot& phi_i) {
    // imex_mod scans all cols if even one is bad. So do that here. The scan is
    // suboptimal, but it's much better to do it like this than have separate
    // kernels just b/c of a safety scan. All scans except the initial guess,
    // which we do indeed handle in a separate kernel, are triggered only very
    // occasionally.
    loop_ki(kv, 1, nvec, [&] (int, int i) {
      for (int k = nlev-1; k >= 0; --k)
        phi_i(k,i) = phi_i(k+1,i) - dphi(k,i);
    });
    return true;
  }

  // Debug routine, not for GPU.
  template <typename R>
  static bool good (const R& a) {
    bool good = true;
    for (int k = 0; k < num_phys_lev; ++k)
      for (int i = 0; i < a.extent_int(1); ++i)
        for (int s = 0; s < packn; ++s) {
          if (i*packn + s >= scaln) continue;
          const auto v = a(k,i)[s];
          if (std::isnan(v) || std::isinf(v)) good = false;
        }
    return good;
  }
};

} // namespace Homme

#endif // HOMMEXX_DIRK_FUNCTOR_IMPL_HPP
