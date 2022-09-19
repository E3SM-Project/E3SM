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
  enum : int { num_work = 12 };

  static_assert(max_num_lev_aligned >= 3,
                "We use wrk(0:2,:) and so need max_num_lev_aligned >= 3");

  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MT = typename TeamPolicy::member_type;

  using Buf1 = ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV_P]>;
  using Buf2 = ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV_P]>;

  using DeparturePoints = ExecViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP][3]>;

  struct Data {
    int nelemd, qsize, limiter_option, cdr_check, hv_q, hv_subcycle_q;
    int geometry_type; // 0: sphere, 1: plane
    Real nu_q, hv_scaling, dp_tol;
    bool independent_time_steps;

    Buf1 buf1[3];
    Buf2 buf2[2];

    DeparturePoints dep_pts;

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
  void reset(const SimulationParams& params);
  int requested_buffer_size() const;
  void init_buffers(const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges();

  void run(const TimeLevel& tl, const Real dt);
  void remap_q(const TimeLevel& tl);

  void calc_trajectory(const int np1, const Real dt);
  void remap_v(const ExecViewUnmanaged<const Scalar*[NUM_TIME_LEVELS][NP][NP][NUM_LEV]>& dp3d,
               const int np1, const ExecViewUnmanaged<const Scalar*[NP][NP][NUM_LEV]>& dp,
               const ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV]>& v);

  void advance_hypervis_scalar(const Real dt);

  int run_trajectory_unit_tests();
  ComposeTransport::TestDepView::HostMirror
  test_trajectory(Real t0, Real t1, const bool independent_time_steps);

  // In test code, the bfb flag says to construct manufactured fields on host to
  // avoid non-bfb-ness in, e.g., trig functions.
  void test_2d(const bool bfb, const int nstep, std::vector<Real>& eval);

  template <int KLIM, typename Fn> KOKKOS_INLINE_FUNCTION
  static void loop_ijk (const KernelVariables& kv, const Fn& h) {
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

  template <typename View> static KOKKOS_INLINE_FUNCTION
  Real* pack2real (const View& v) { return &(*v.data())[0]; }
  template <typename View> static KOKKOS_INLINE_FUNCTION
  const Real* cpack2real (const View& v) { return &(*v.data())[0]; }
};

} // namespace Homme

#endif // HOMMEXX_COMPOSE_TRANSPORT_IMPL_HPP
#endif // HOMME_ENABLE_COMPOSE
