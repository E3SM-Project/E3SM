/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#ifndef HOMMEXX_COMPOSE_TRANSPORT_IMPL_HPP
#define HOMMEXX_COMPOSE_TRANSPORT_IMPL_HPP

#include "ComposeTransport.hpp"

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

struct ComposeTransportImpl {
  enum : int { np = NP };
  enum : int { packn = VECTOR_SIZE };
  enum : int { np2 = NP*NP };
  enum : int { num_lev_pack = NUM_LEV };
  enum : int { max_num_lev_pack = NUM_LEV_P };
  enum : int { max_num_lev_aligned = max_num_lev_pack*packn };
  enum : int { num_phys_lev = NUM_PHYSICAL_LEV };

  static_assert(max_num_lev_aligned >= 3,
                "We use wrk(0:2,:) and so need max_num_lev_aligned >= 3");

  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MT = typename TeamPolicy::member_type;

  // For the enhanced trajectory, we need one extra level beyond the usual.
  using Buf1Alloc = ExecViewUnmanaged<Scalar*[NP][NP][
    ((NUM_PHYSICAL_LEV + 2) + VECTOR_SIZE - 1) / VECTOR_SIZE]>;
  using Buf1o = ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV_P]>;
  using Buf1e = Buf1Alloc;

  using Buf2 = ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV_P]>;

  using DeparturePoints = ExecViewManaged<Real*****>;

  typedef ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]> SNlev;
  typedef ExecViewUnmanaged<Real[NP][NP][NUM_LEV*VECTOR_SIZE]> RNlev;
  typedef ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV_P]> SNlevp;
  typedef ExecViewUnmanaged<Real[NP][NP][NUM_LEV_P*VECTOR_SIZE]> RNlevp;
  typedef ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]> S2Nlev;
  typedef ExecViewUnmanaged<Real[2][NP][NP][NUM_LEV*VECTOR_SIZE]> R2Nlev;
  typedef ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV_P]> S2Nlevp;
  typedef typename ViewConst<SNlev >::type CSNlev;
  typedef typename ViewConst<RNlev >::type CRNlev;
  typedef typename ViewConst<SNlevp>::type CSNlevp;
  typedef typename ViewConst<RNlevp>::type CRNlevp;
  typedef typename ViewConst<S2Nlev>::type CS2Nlev;
  typedef typename ViewConst<R2Nlev>::type CR2Nlev;

  using  DpSlot = ExecViewUnmanaged<      Scalar**   [NP][NP][NUM_LEV]>;
  using   VSlot = ExecViewUnmanaged<      Scalar**[2][NP][NP][NUM_LEV]>;
  using CDpSlot = ExecViewUnmanaged<const Scalar**   [NP][NP][NUM_LEV]>;
  using  CVSlot = ExecViewUnmanaged<const Scalar**[2][NP][NP][NUM_LEV]>;
  struct VelocityRecord;

  struct Data {
    int nelemd, qsize, limiter_option, cdr_check, hv_q, hv_subcycle_q;
    int geometry_type; // 0: sphere, 1: plane
    int trajectory_nsubstep; // 0: original alg, >= 1: enhanced
    Real nu_q, hv_scaling, dp_tol, deta_tol;
    bool independent_time_steps;

    // buf1o and buf1e point to the same memory, sized to the larger of the
    // two. They are used in different parts of the code.
    static constexpr int n_buf1 = 4, n_buf2 = 4;
    Buf1o buf1o[n_buf1];
    Buf1e buf1e[n_buf1];
    Buf2 buf2[n_buf2];

    ExecView<Scalar[NUM_LEV]> hydetai; // diff(etai)
    ExecView<Real[NUM_INTERFACE_LEV]> hydetam_ref;

    DeparturePoints dep_pts, vnode, vdep; // (ie,lev,i,j,d)

    std::shared_ptr<VelocityRecord> vrec;

    Data ()
      : nelemd(-1), qsize(-1), limiter_option(9), cdr_check(0), hv_q(0),
        hv_subcycle_q(0), geometry_type(0), nu_q(0), hv_scaling(0), dp_tol(-1),
        independent_time_steps(false)
    {}
  };

  HybridVCoord m_hvcoord;
  Elements m_elements;
  ElementsState m_state;
  ElementsDerivedState m_derived;
  ElementsGeometry m_geometry;
  Tracers m_tracers;
  SphereOperators m_sphere_ops;
  int nslot;
  Data m_data;

  TeamPolicy m_tp_ne, m_tp_ne_qsize, m_tp_ne_hv_q;
  TeamUtils<ExecSpace> m_tu_ne, m_tu_ne_qsize, m_tu_ne_hv_q;

  std::shared_ptr<BoundaryExchange>
    m_qdp_dss_be[Q_NUM_TIME_LEVELS], m_v_dss_be[2], m_hv_dss_be[2];

  ComposeTransportImpl();
  ComposeTransportImpl(const int num_elems);

  void setup();

  KOKKOS_INLINE_FUNCTION
  size_t shmem_size (const int team_size) const {
    return KernelVariables::shmem_size(team_size);
  }

  void set_dp_tol();
  void setup_enhanced_trajectory();
  void reset(const SimulationParams& params);
  int requested_buffer_size() const;
  void init_buffers(const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges();

  void run(const TimeLevel& tl, const Real dt);
  void remap_q(const TimeLevel& tl);

  void calc_trajectory(const int np1, const Real dt);
  void calc_enhanced_trajectory(const int np1, const Real dt);
  void remap_v(const ExecViewUnmanaged<const Scalar*[NUM_TIME_LEVELS][NP][NP][NUM_LEV]>& dp3d,
               const int np1, const ExecViewUnmanaged<const Scalar*[NP][NP][NUM_LEV]>& dp,
               const ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV]>& v);

  void advance_hypervis_scalar(const Real dt);

  int run_trajectory_unit_tests();
  int run_enhanced_trajectory_unit_tests();
  ComposeTransport::TestDepView::HostMirror
  test_trajectory(Real t0, Real t1, const bool independent_time_steps);

  // In test code, the bfb flag says to construct manufactured fields on host to
  // avoid non-bfb-ness in, e.g., trig functions.
  void test_2d(const bool bfb, const int nstep, std::vector<Real>& eval);

  template <typename Fn> KOKKOS_INLINE_FUNCTION
  static void loop_ijk (const int KLIM, const KernelVariables& kv, const Fn& h) {
    using Kokkos::parallel_for;
    if (OnGpu<ExecSpace>::value) {
      const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
      const auto tvr = Kokkos::ThreadVectorRange(kv.team, KLIM);
      const auto f = [&] (const int idx) {
        const int i = idx / NP, j = idx % NP;
        const auto g = [&] (const int k) { h(i,j,k); };
        parallel_for(tvr, g);
      };
      parallel_for(ttr, f);
    } else if (kv.team.team_size() == 1) {
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < KLIM; ++k)
            h(i,j,k);
    } else {
      const auto tr = Kokkos::TeamThreadRange(kv.team, KLIM);
      const auto f = [&] (const int k) {
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j)
            h(i,j,k);
      };
      parallel_for(tr, f);
    }
  }

  template <int KLIM, typename Fn> KOKKOS_INLINE_FUNCTION
  static void loop_ijk (const KernelVariables& kv, const Fn& h) {
    loop_ijk(KLIM, kv, h);
  }

  template <typename Fn> KOKKOS_INLINE_FUNCTION
  static void loop_ij (const KernelVariables& kv, const Fn& h) {
    if (OnGpu<ExecSpace>::value) {
      const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
      const auto f = [&] (const int idx) {
        const int i = idx / NP, j = idx % NP;
        const auto g = [&] () { h(i,j); };
        Kokkos::single(Kokkos::PerThread(kv.team), g);
      };
      Kokkos::parallel_for(ttr, f);
    } else if (kv.team.team_size() == 1) {
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          h(i,j);
    } else {
      const auto f = [&] () {
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j)
            h(i,j);
      };
      Kokkos::single(Kokkos::PerTeam(kv.team), f);
    }
  }

  template <typename Fn>
  void loop_host_ie_plev_ij (const Fn& f) const {
    for (int ie = 0; ie < m_data.nelemd; ++ie)
      for (int lev = 0; lev < num_phys_lev; ++lev)
        for (int i = 0; i < np; ++i)
          for (int j = 0; j < np; ++j)
            f(ie, lev, i, j);
  }

  template <int nlev> KOKKOS_INLINE_FUNCTION
  static void idx_ie_nlev_ij (const int idx, int& ie, int& lev, int& i, int& j) {
    ie = idx / (nlev*np*np);
    lev = (idx / (np*np)) % nlev;
    i = (idx / np) % np;
    j = idx % np;
  }

  KOKKOS_INLINE_FUNCTION
  static void idx_ie_physlev_ij (const int idx, int& ie, int& lev, int& i, int& j) {
    idx_ie_nlev_ij<num_phys_lev>(idx, ie, lev, i, j);
  }

  template <typename Fn>
  void launch_ie_physlev_ij (Fn& f) const {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, m_data.nelemd*np*np*num_phys_lev), f);
  }

  KOKKOS_INLINE_FUNCTION
  static void idx_ie_packlev_ij (const int idx, int& ie, int& lev, int& i, int& j) {
    idx_ie_nlev_ij<num_lev_pack>(idx, ie, lev, i, j);
  }

  template <typename Fn>
  void launch_ie_packlev_ij (Fn& f) const {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, m_data.nelemd*np*np*num_lev_pack), f);
  }

  template <int nlev> KOKKOS_INLINE_FUNCTION
  static void idx_ie_ij_nlev (const int idx, int& ie, int& i, int& j, int& lev) {
    ie = idx / (np*np*nlev);
    i = (idx / (np*nlev)) % np;
    j = (idx / nlev) % np;
    lev = idx % nlev;
  }

  template <int nlev, typename Fn>
  void launch_ie_ij_nlev (Fn& f) const {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, m_data.nelemd*np*np*nlev), f);
  }

  template <int nlev> KOKKOS_INLINE_FUNCTION
  static void idx_ie_q_ij_nlev (const int qsize, const int idx,
                                int& ie, int& q, int& i, int& j, int& lev) {
    ie = idx / (qsize*np*np*nlev);
    q = (idx / (np*np*nlev)) % qsize;
    i = (idx / (np*nlev)) % np;
    j = (idx / nlev) % np;
    lev = idx % nlev;
  }

  template <int nlev, typename Fn>
  void launch_ie_q_ij_nlev (const int qsize, Fn& f) const {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, m_data.nelemd*qsize*np*np*nlev), f);
  }

  template <typename V>
  static decltype(Kokkos::create_mirror_view(V())) cmvdc (const V& v) {
    const auto h = Kokkos::create_mirror_view(v);
    deep_copy(h, v);
    return h;
  }

  static KOKKOS_INLINE_FUNCTION
  Real* pack2real (Scalar* pack) { return &(*pack)[0]; }
  static KOKKOS_INLINE_FUNCTION
  const Real* pack2real (const Scalar* pack) { return &(*pack)[0]; }
  template <typename View> static KOKKOS_INLINE_FUNCTION
  Real* pack2real (const View& v) { return pack2real(v.data()); }
  template <typename View> static KOKKOS_INLINE_FUNCTION
  const Real* cpack2real (const View& v) { return pack2real(v.data()); }

  KOKKOS_FUNCTION
  static void ugradv_sphere (
    const SphereOperators& sphere_ops, const KernelVariables& kv,
    const typename ViewConst<ExecViewUnmanaged<Real[2][3][NP][NP]> >::type& vec_sphere2cart,
    // velocity, latlon
    const typename ViewConst<ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]> >::type& u,
    const typename ViewConst<ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]> >::type& v,
    const ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>& v_cart,
    const ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]>& ugradv_cart,
    // [u dot grad] v, latlon
    const ExecViewUnmanaged<Scalar[2][NP][NP][NUM_LEV]>& ugradv)
  {
    for (int d_cart = 0; d_cart < 3; ++d_cart) {
      const auto f1 = [&] (const int i, const int j, const int k) {
        v_cart(i,j,k) = (vec_sphere2cart(0,d_cart,i,j) * v(0,i,j,k) +
                         vec_sphere2cart(1,d_cart,i,j) * v(1,i,j,k));      
      };
      loop_ijk<NUM_LEV>(kv, f1);
      kv.team_barrier();

      sphere_ops.gradient_sphere<NUM_LEV>(kv, v_cart, ugradv_cart);

      const auto f2 = [&] (const int i, const int j, const int k) {
        if (d_cart == 0) ugradv(0,i,j,k) = ugradv(1,i,j,k) = 0;
        for (int d_latlon = 0; d_latlon < 2; ++d_latlon)
          ugradv(d_latlon,i,j,k) +=
            vec_sphere2cart(d_latlon,d_cart,i,j)*
            (u(0,i,j,k) * ugradv_cart(0,i,j,k) + u(1,i,j,k) * ugradv_cart(1,i,j,k));
      };
      loop_ijk<NUM_LEV>(kv, f2);
    }
  }

  // Form a 3rd-degree Lagrange polynomial over (x(k-1:k+1), y(k-1:k+1)) and set
  // yi(k) to its derivative at x(k). yps(:,:,0) is not written.
  template <typename Real>
  KOKKOS_FUNCTION static Real approx_derivative (
    const Real& xkm1, const Real& xk, const Real& xkp1,
    const Real& ykm1, const Real& yk, const Real& ykp1)
  {
    return (ykm1*((         1 /(xkm1 - xk  ))*((xk - xkp1)/(xkm1 - xkp1))) +
            yk  *((         1 /(xk   - xkm1))*((xk - xkp1)/(xk   - xkp1)) +
                  ((xk - xkm1)/(xk   - xkm1))*(         1 /(xk   - xkp1))) +
            ykp1*(((xk - xkm1)/(xkp1 - xkm1))*(         1 /(xkp1 - xk  ))));
  }

  KOKKOS_INLINE_FUNCTION static void approx_derivative (
    const KernelVariables& kv, const CSNlevp& xs, const CSNlevp& ys,
    const SNlev& yps) // yps(:,:,0) is undefined
  {
    CRNlevp x(cpack2real(xs));
    CRNlevp y(cpack2real(ys));
    RNlev yp(pack2real(yps));
    const auto f = [&] (const int i, const int j, const int k) {
      if (k == 0) return;
      yp(i,j,k) = approx_derivative(x(i,j,k-1), x(i,j,k), x(i,j,k+1),
                                    y(i,j,k-1), y(i,j,k), y(i,j,k+1));
    };
    loop_ijk<num_phys_lev>(kv, f);
  }

  template <typename HyBiPackT, typename DivDpScalT,
            typename EddPackT, typename EddScalT>
  KOKKOS_FUNCTION static void calc_eta_dot_dpdn (
    const KernelVariables& kv,
    const HyBiPackT& hybrid_bi, // const Scalar[NUM_LEV_P]
    // divergence_sphere of (v dp) at midpoints, scalar
    const DivDpScalT& divdps,
    // eta_dot_dpdn at interfaces, pack and scalar views of same data
    const EddPackT& edd, const EddScalT& edds)
  {
    const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
    const auto tvr = Kokkos::ThreadVectorRange(kv.team, NUM_LEV);
    const auto f1 = [&] (const int idx) {
      const int i = idx / NP, j = idx % NP;
      const auto r = [&] (const int k, Real& dps, const bool final) {
        assert(k != 0 || dps == 0);
        if (final) edds(i,j,k) = dps;
        dps += divdps(i,j,k);
      };
      Dispatch<>::parallel_scan(kv.team, num_phys_lev, r);
    };
    Kokkos::parallel_for(ttr, f1);
    kv.team_barrier();
    const auto f2 = [&] (const int idx) {
      const int i = idx / NP, j = idx % NP;
      const int kend = num_phys_lev - 1;
      const Real dps = edds(i,j,kend) + divdps(i,j,kend);
      assert(hybrid_bi(0)[0] == 0);
      const auto s = [&] (const int kp) {
        edd(i,j,kp) = hybrid_bi(kp)*dps - edd(i,j,kp);
        if (kp == 0) edd(i,j,kp)[0] = 0;
      };
      Kokkos::parallel_for(tvr, s);
    };
    Kokkos::parallel_for(ttr, f2);
    kv.team_barrier();
    const int bottom = num_phys_lev;
    const auto f3 = [&] (const int idx) {
      const int i = idx / NP, j = idx % NP;
      Kokkos::single(Kokkos::PerThread(kv.team), [&] () { edds(i,j,bottom) = 0; });
    };
    Kokkos::parallel_for(ttr, f3);
  }
};

} // namespace Homme

#endif // HOMMEXX_COMPOSE_TRANSPORT_IMPL_HPP
#endif // HOMME_ENABLE_COMPOSE
