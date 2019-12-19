/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "ColumnOps.hpp"
#include "ElementsDerivedState.hpp"
#include "KernelVariables.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/TestUtils.hpp"

#include <limits>
#include <random>

namespace Homme {

void ElementsDerivedState::init(const int num_elems) {
  // Sanity check
  assert (num_elems>0);

  m_num_elems = num_elems;

  m_omega_p = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Omega P", m_num_elems);

  m_vn0 = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Derived Lateral Velocities", m_num_elems);

  m_eta_dot_dpdn = ExecViewManaged<Scalar * [NP][NP][NUM_LEV_P]>("eta_dot_dpdn", m_num_elems);

  m_dp                = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dp", m_num_elems);
  m_divdp             = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_divdp", m_num_elems);
  m_divdp_proj        = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_divdp_proj", m_num_elems);
  m_dpdiss_biharmonic = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dpdiss_biharmonic", m_num_elems);
  m_dpdiss_ave        = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dpdiss_ave", m_num_elems);
}

void ElementsDerivedState::randomize(const int seed, const Real dp3d_min) {
  // Check derived state was inited
  assert (m_num_elems>0);

  // Sanity check
  assert (dp3d_min>0);

  // arbitrary minimum value to generate
  constexpr const Real min_value = 0.015625;
  std::mt19937_64 engine(seed);
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(m_omega_p, engine, random_dist);
  genRandArray(m_vn0,     engine, random_dist);

  // Generate eta_dot_dpdn so that it is << dp3d
  genRandArray(m_eta_dot_dpdn, engine, std::uniform_real_distribution<Real>(0.01*dp3d_min,0.1*dp3d_min));

  // During remap, in the case of rsplit=0, the target thickness is computed as
  //   dp_star(k) = dp(k) + dt*(eta_dot_dpdn(k+1)-eta_dot_dpdn(k))
  // Such dp_start should verify sum(dp_star) = sum(dp). One way to enforce this
  // is to modify one entry of eta_dot_dpdn, so that the two are "close enough"
  // to not raise errors during unit tests (inside unit tests we check the constraint
  // with a somewhat loose tolerance).
  auto eta_dot_dpdn = m_eta_dot_dpdn;
  auto policy = Homme::get_default_team_policy<ExecSpace>(m_num_elems);
  const int nteams = Homme::get_num_concurrent_teams(policy);
  ExecViewManaged<Scalar *[NP][NP][NUM_LEV]> delta_eta("",nteams);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const TeamMember& team) {
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto eta_dot_dpdn_ij = Homme::subview(eta_dot_dpdn,kv.ie,igp,jgp);
      auto delta_eta_ij    = Homme::subview(delta_eta,kv.team_idx,igp,jgp);

      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        constexpr int LAST_PACK     = ColInfo<NUM_INTERFACE_LEV>::LastPack;
        constexpr int LAST_PACK_LEN = ColInfo<NUM_INTERFACE_LEV>::LastPackLen;
        constexpr int LAST_PACK_END = ColInfo<NUM_INTERFACE_LEV>::LastPackEnd;

        // Last entry of eta_dot_dpdn should be 0
        eta_dot_dpdn_ij(LAST_PACK)[LAST_PACK_END] = 0.0;

        // If VECTOR_SIZE does not divide NUM_INTERFACE_LEV, zero out the
        // garbage at the end
        if (LAST_PACK_LEN!=VECTOR_SIZE) {
          for (int i=LAST_PACK_LEN; i<VECTOR_SIZE; ++i) {
            eta_dot_dpdn_ij(LAST_PACK)[i] = 0.0;
          }
        }
      });

      // Compute delta_eta over the column
      ColumnOps::compute_midpoint_delta(kv,eta_dot_dpdn_ij,delta_eta_ij);

      Real sum = 0;
      // Integrate delta_eta over the column, and subtract that value
      // to the first entry of eta_dot_dpdn
      ColumnOps::column_reduction<NUM_PHYSICAL_LEV>(kv.team,delta_eta_ij,sum);

      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        eta_dot_dpdn_ij(0)[0] -= sum;
      });
    });
  });
}

} // namespace Homme
