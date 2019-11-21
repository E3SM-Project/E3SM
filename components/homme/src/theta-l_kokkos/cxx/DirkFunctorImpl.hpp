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
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "SimulationParams.hpp"
#include "PhysicalConstants.hpp"
#include "profiling.hpp"
#include "ErrorDefs.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/ViewUtils.hpp"
#include "utilities/scream_tridiag.hpp"

#include <assert.h>

namespace Homme {

struct DirkFunctorImpl {
  enum : int { packn = VECTOR_SIZE };
  enum : int { scaln = NP*NP };
  enum : int { npack = (scaln + packn - 1)/packn };
  enum : int { max_num_lev_pack = NUM_LEV_P };
  enum : int { num_lev_aligned = max_num_lev_pack*packn };
  enum : int { num_phys_lev = NUM_PHYSICAL_LEV };
  enum : int { num_work = 8 };

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
  TeamPolicy m_policy;

  KOKKOS_INLINE_FUNCTION
  size_t shmem_size (const int team_size) const {
    return KernelVariables::shmem_size(team_size);
  }

  DirkFunctorImpl (const int& nelem)
    : m_policy(1,1,1) // throwaway settings
  {
    init(nelem);
  }

  void init (const int& nelem) {
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
      m_policy = TeamPolicy(nelem, 1, 1);
    }
    const int nteam = get_num_concurrent_teams(m_policy);
    m_work = Work("DirkFunctorImpl::m_work", nteam);
    m_ls = LinearSystem("DirkFunctorImpl::m_ls", nteam);
  }

  void run () {
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

    const int n = npack;
    const auto pv = Kokkos::ThreadVectorRange(kv.team, n);

    const auto f1 = [&] (const int) {
      const auto k0 = [&] (const int i) { gwh_i(0,i) = 0; };
      parallel_for(pv, k0);
    };
    parallel_for(Kokkos::TeamThreadRange(kv.team, 1), f1);

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
    parallel_for(Kokkos::TeamThreadRange(kv.team, nlev-1), f2);
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
    using Kokkos::subview;
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
};

} // Namespace Homme

#endif // HOMMEXX_DIRK_FUNCTOR_IMPL_HPP
