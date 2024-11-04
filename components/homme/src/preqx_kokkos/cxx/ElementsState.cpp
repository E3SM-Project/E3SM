/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "ElementsState.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"
#include "HybridVCoord.hpp"

#include <limits>
#include <random>
#include <assert.h>

namespace Homme {

void ElementsState::init(const int num_elems) {
  // Sanity check
  assert (num_elems>0);

  m_num_elems = num_elems;

  m_v    = ExecViewManaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV]>("Horizontal Velocity", m_num_elems);
  m_t    = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV]>("Temperature", m_num_elems);
  m_dp3d = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV]>("DP3D", m_num_elems);

  m_ps_v = ExecViewManaged<Real * [NUM_TIME_LEVELS][NP][NP]>("PS_V", m_num_elems);
}

void ElementsState::randomize(const int seed) {
  randomize(seed,1.0);
}

void ElementsState::randomize(const int seed, const Real max_pressure) {
  randomize(seed,max_pressure,max_pressure/100, 0.0);
}

void ElementsState::randomize(const int seed, const Real max_pressure, const Real ps0, const Real hyai0) {
  // Check state was inited
  assert (m_num_elems>0);

  // Check data makes sense
  assert (max_pressure>ps0);
  assert (ps0>0);
  assert (hyai0>=0);

  // arbitrary minimum value to generate and minimum determinant allowed
  constexpr const Real min_value = 0.015625;
  std::mt19937_64 engine(seed);
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(m_v, engine, random_dist);
  genRandArray(m_t, engine, random_dist);

  // This ensures the pressure in a single column is monotonically increasing
  // and has fixed upper and lower values
  const auto make_pressure_partition = [=](
      ExecViewUnmanaged<Scalar[NUM_LEV]> pt_dp) {

    auto h_pt_dp = Kokkos::create_mirror_view(pt_dp);
    Kokkos::deep_copy(h_pt_dp,pt_dp);
    Real* data     = reinterpret_cast<Real*>(h_pt_dp.data());
    Real* data_end = data + NUM_PHYSICAL_LEV;

    Real p[NUM_INTERFACE_LEV];
    Real* p_start = &p[0];
    Real* p_end   = p_start+NUM_INTERFACE_LEV;

    for (int i=0; i<NUM_PHYSICAL_LEV; ++i) {
      p[i+1] = data[i];
    }
    p[0] = ps0;
    p[NUM_INTERFACE_LEV-1] = max_pressure;

    // Put in monotonic order
    std::sort(p_start, p_end);

    // Check for no repetitions
    if (std::unique(p_start,p_end)!=p_end) {
      return false;
    }

    // Compute dp from p (we assume p(last interface)=max_pressure)
    for (int i=0; i<NUM_PHYSICAL_LEV; ++i) {
      data[i] = p[i+1]-p[i];
    }

    // Check that dp>=dp_min
    const Real min_dp = std::numeric_limits<Real>::epsilon()*1000;
    for (auto it=data; it!=data_end; ++it) {
      if (*it < min_dp) {
        return false;
      }
    }

    // Fill remainder of last vector pack with quiet nan's
    Real* alloc_end = data+NUM_LEV*VECTOR_SIZE;
    for (auto it=data_end; it!=alloc_end; ++it) {
      *it = std::numeric_limits<Real>::quiet_NaN();
    }

    Kokkos::deep_copy(pt_dp,h_pt_dp);

    return true;
  };

  std::uniform_real_distribution<Real> pressure_pdf(ps0, max_pressure);

  for (int ie = 0; ie < m_num_elems; ++ie) {
    // Because this constraint is difficult to satisfy for all of the tensors,
    // incrementally generate the view
    for (int igp = 0; igp < NP; ++igp) {
      for (int jgp = 0; jgp < NP; ++jgp) {
        for (int tl = 0; tl < NUM_TIME_LEVELS; ++tl) {
          ExecViewUnmanaged<Scalar[NUM_LEV]> pt_dp3d =
              Homme::subview(m_dp3d, ie, tl, igp, jgp);
          do {
            genRandArray(pt_dp3d, engine, pressure_pdf);
          } while (make_pressure_partition(pt_dp3d)==false);
        }
      }
    }
  }

  // Generate ps_v so that it is equal to sum(dp3d).
  HybridVCoord hvcoord;
  hvcoord.ps0 = ps0;
  hvcoord.hybrid_ai0 = hyai0;
  hvcoord.m_inited = true;
  auto dp = m_dp3d;
  auto ps = m_ps_v;
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(m_num_elems*NUM_TIME_LEVELS),
                       KOKKOS_LAMBDA(const TeamMember team){
    KernelVariables kv(team);
    const int ie = kv.ie / NUM_TIME_LEVELS;
    const int tl = kv.ie % NUM_TIME_LEVELS;
    hvcoord.compute_ps_ref_from_dp(kv,Homme::subview(dp,ie,tl),
                                      Homme::subview(ps,ie,tl));
  });
}

void ElementsState::pull_from_f90_pointers (CF90Ptr& state_v,    CF90Ptr& state_t,
                                            CF90Ptr& state_dp3d, CF90Ptr& state_ps_v) {
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_t_f90    (state_t,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_dp3d_f90 (state_dp3d,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]> state_v_f90    (state_v,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS]                     [NP][NP]> ps_v_f90       (state_ps_v,m_num_elems);

  sync_to_device(state_t_f90,    m_t);
  sync_to_device(state_dp3d_f90, m_dp3d);
  sync_to_device(state_v_f90,    m_v);

  // F90 ptrs to arrays (np,np,num_time_levels,nelemd) can be stuffed directly in an unmanaged view
  // with scalar Real*[NUM_TIME_LEVELS][NP][NP] (with runtime dimension nelemd)

  auto ps_v_host = Kokkos::create_mirror_view(m_ps_v);
  Kokkos::deep_copy(ps_v_host,ps_v_f90);
  Kokkos::deep_copy(m_ps_v,ps_v_host);
}

void ElementsState::push_to_f90_pointers (F90Ptr& state_v, F90Ptr& state_t, F90Ptr& state_dp3d) const {
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_t_f90    (state_t,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV]   [NP][NP]> state_dp3d_f90 (state_dp3d,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]> state_v_f90    (state_v,m_num_elems);

  sync_to_host(m_t,    state_t_f90);
  sync_to_host(m_dp3d, state_dp3d_f90);
  sync_to_host(m_v,    state_v_f90);
}

void check_print_abort_on_bad_elems (const std::string& label, const int time_level,
                                     const int error_code) {
  // Not implemented.
  return;
}

HashType ElementsState::hash (const int) const { return 0; }

} // namespace Homme
