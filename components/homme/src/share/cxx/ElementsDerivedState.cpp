/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "ElementsDerivedState.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"

#include <limits>
#include <random>

namespace Homme {

void ElementsDerivedState::init(const int num_elems) {
  m_num_elems = num_elems;

  m_omega_p = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Omega P", m_num_elems);

  m_vn0 = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Derived Lateral Velocities", m_num_elems);

  m_eta_dot_dpdn = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("eta_dot_dpdn", m_num_elems);

  m_dp                = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dp", m_num_elems);
  m_divdp             = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_divdp", m_num_elems);
  m_divdp_proj        = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_divdp_proj", m_num_elems);
  m_dpdiss_biharmonic = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dpdiss_biharmonic", m_num_elems);
  m_dpdiss_ave        = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dpdiss_ave", m_num_elems);
}

void ElementsDerivedState::random_init(int num_elems, Real dp3d_min) {
  init(num_elems);

  // arbitrary minimum value to generate
  constexpr const Real min_value = 0.015625;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(m_omega_p, engine, random_dist);
  genRandArray(m_vn0,     engine, random_dist);

  // Generate eta_dot_dpdn so that it is << dp3d
  genRandArray(m_eta_dot_dpdn, engine, std::uniform_real_distribution<Real>(0.01*dp3d_min,0.1*dp3d_min));
}

} // namespace Homme
