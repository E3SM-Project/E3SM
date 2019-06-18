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

void StateStorage::init_storage(const int num_elems) {

  m_v         = ExecViewManaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV  ]>("Horizontal velocity", num_elems);
  m_w_i       = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV_P]>("Vertical velocity at interfaces", num_elems);
  m_vtheta_dp = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV  ]>("Virtual potential temperature", num_elems);
  m_phinh_i   = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV_P]>("Geopotential at interfaces", num_elems);
  m_dp3d      = ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV  ]>("Delta p at levels", num_elems);

  m_ps_v = ExecViewManaged<Real * [NUM_TIME_LEVELS][NP][NP]>("PS_V", num_elems);
}

void StateStorage::copy_state(const StateStorage& src) {

  Kokkos::deep_copy(m_v        , src.m_v         );
  Kokkos::deep_copy(m_w_i      , src.m_w_i       );
  Kokkos::deep_copy(m_vtheta_dp, src.m_vtheta_dp );
  Kokkos::deep_copy(m_phinh_i  , src.m_phinh_i   );
  Kokkos::deep_copy(m_dp3d     , src.m_dp3d      );
  Kokkos::deep_copy(m_ps_v     , src.m_ps_v      );
}

void ElementsState::init(const int num_elems) {
  m_num_elems = num_elems;
  StateStorage::init_storage(num_elems);
}

//test for tensor hv is needed
void ElementsState::random_init(const int num_elems, const int seed, const Real max_pressure) {
  HybridVCoord hv;
  hv.random_init(seed);
  random_init(num_elems,seed,max_pressure,hv);
}

void ElementsState::random_init(const int num_elems, const int seed, 
                                const Real max_pressure, const HybridVCoord& hvcoord) {
  // arbitrary minimum value to generate and minimum determinant allowed
  constexpr const Real min_value = 0.015625;

  // We may re-init elements in a unit test. If so, don't reallocate, since it could mess up
  // the Context structure.
  if (m_num_elems==0) {
    init(num_elems);
  }
  std::mt19937_64 engine(seed);
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(m_v,         engine, random_dist);
  genRandArray(m_w_i,       engine, random_dist);
  genRandArray(m_vtheta_dp, engine, random_dist);
  genRandArray(m_phinh_i,   engine, random_dist);

  // Generate ps_v so that it is >> ps0.
  // Note: make sure you init hvcoord before calling this method!
  genRandArray(m_ps_v, engine, std::uniform_real_distribution<Real>(100*hvcoord.ps0,1000*hvcoord.ps0));

  // This ensures the pressure in a single column is monotonically increasing
  // and has fixed upper and lower values
  const auto make_pressure_partition = [=](
      HostViewUnmanaged<Scalar[NUM_LEV]> pt_pressure) {
    // Put in monotonic order
    std::sort(
        reinterpret_cast<Real *>(pt_pressure.data()),
        reinterpret_cast<Real *>(pt_pressure.data() + pt_pressure.size()));
    // Ensure none of the values are repeated
    for (int level = NUM_PHYSICAL_LEV - 1; level > 0; --level) {
      const int prev_ilev = (level - 1) / VECTOR_SIZE;
      const int prev_vlev = (level - 1) % VECTOR_SIZE;
      const int cur_ilev = level / VECTOR_SIZE;
      const int cur_vlev = level % VECTOR_SIZE;
      // Need to try again if these are the same or if the thickness is too
      // small
      if (pt_pressure(cur_ilev)[cur_vlev] <=
          pt_pressure(prev_ilev)[prev_vlev] +
              min_value * std::numeric_limits<Real>::epsilon()) {
        return false;
      }
    }
    // We know the minimum thickness of a layer is min_value * epsilon
    // (due to floating point), so set the bottom layer thickness to that,
    // and subtract that from the top layer
    // This ensures that the total sum is max_pressure
    pt_pressure(0)[0] = min_value * std::numeric_limits<Real>::epsilon();
    const int top_ilev = (NUM_PHYSICAL_LEV - 1) / VECTOR_SIZE;
    const int top_vlev = (NUM_PHYSICAL_LEV - 1) % VECTOR_SIZE;
    // Note that this may not actually change the top level pressure
    // This is okay, because we only need to approximately sum to max_pressure
    pt_pressure(top_ilev)[top_vlev] = max_pressure - pt_pressure(0)[0];
    for (int e_vlev = top_vlev + 1; e_vlev < VECTOR_SIZE; ++e_vlev) {
      pt_pressure(top_ilev)[e_vlev] = std::numeric_limits<Real>::quiet_NaN();
    }
    // Now compute the interval thicknesses
    for (int level = NUM_PHYSICAL_LEV - 1; level > 0; --level) {
      const int prev_ilev = (level - 1) / VECTOR_SIZE;
      const int prev_vlev = (level - 1) % VECTOR_SIZE;
      const int cur_ilev = level / VECTOR_SIZE;
      const int cur_vlev = level % VECTOR_SIZE;
      pt_pressure(cur_ilev)[cur_vlev] -= pt_pressure(prev_ilev)[prev_vlev];
    }
    return true;
  };

  std::uniform_real_distribution<Real> pressure_pdf(min_value, max_pressure);

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
}

void ElementsState::save_state ()
{
  if (m_state0.m_v.extent_int(0)==0) {
    m_state0.init_storage(m_num_elems);
  }
  m_state0.copy_state(*this);
}

void ElementsState::pull_from_f90_pointers (CF90Ptr& state_v,         CF90Ptr& state_w_i,
                                            CF90Ptr& state_vtheta_dp, CF90Ptr& state_phinh_i,
                                            CF90Ptr& state_dp3d,      CF90Ptr& state_ps_v) {
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ][2][NP][NP]> state_v_f90         (state_v,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_w_i_f90       (state_w_i,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_vtheta_dp_f90 (state_vtheta_dp,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_phinh_i_f90   (state_phinh_i,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_dp3d_f90      (state_dp3d,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_TIME_LEVELS]                      [NP][NP]> ps_v_f90            (state_ps_v,m_num_elems);

  sync_to_device(state_v_f90,         m_v);
  sync_to_device(state_w_i_f90,       m_w_i);
  sync_to_device(state_vtheta_dp_f90, m_vtheta_dp);
  sync_to_device(state_phinh_i_f90,   m_phinh_i);
  sync_to_device(state_dp3d_f90,      m_dp3d);

  // F90 ptrs to arrays (np,np,num_time_levels,nelemd) can be stuffed directly in an unmanaged view
  // with scalar Real*[NUM_TIME_LEVELS][NP][NP] (with runtime dimension nelemd)

  auto ps_v_host = Kokkos::create_mirror_view(m_ps_v);
  Kokkos::deep_copy(ps_v_host,ps_v_f90);
  Kokkos::deep_copy(m_ps_v,ps_v_host);
}

void ElementsState::push_to_f90_pointers (F90Ptr& state_v, F90Ptr& state_w_i, F90Ptr& state_vtheta_dp,
                                          F90Ptr& state_phinh_i, F90Ptr& state_dp3d) const {
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ][2][NP][NP]> state_v_f90         (state_v,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_w_i_f90       (state_w_i,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_vtheta_dp_f90 (state_vtheta_dp,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_INTERFACE_LEV]   [NP][NP]> state_phinh_i_f90   (state_phinh_i,m_num_elems);
  HostViewUnmanaged<Real *[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV ]   [NP][NP]> state_dp3d_f90      (state_dp3d,m_num_elems);

  sync_to_host(m_v,         state_v_f90);
  sync_to_host(m_w_i,       state_w_i_f90);
  sync_to_host(m_vtheta_dp, state_vtheta_dp_f90);
  sync_to_host(m_phinh_i,   state_phinh_i_f90);
  sync_to_host(m_dp3d,      state_dp3d_f90);
}

} // namespace Homme
