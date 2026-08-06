#include "ElementsForcing.hpp"

#include "utilities/TestUtils.hpp"
#include <random>

namespace Homme {

void ElementsForcing::init (const int num_elems) {
  // Sanity check
  assert (num_elems>0);

  m_num_elems = num_elems;

  m_fm      = ExecViewManaged<Scalar * [3][NP][NP][NUM_LEV]  >("F_Momentum",    m_num_elems);
  m_fvtheta = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("F_VirtualPotentialTemperature", m_num_elems);
  m_ft      = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("F_Temperature", m_num_elems);
  m_fphi    = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV_P]>("F_Phi",         m_num_elems);
}

void ElementsForcing::randomize (const int seed, const Real min_f, const Real max_f) {
  // Check forcing was inited
  assert (m_num_elems>0);
  assert (max_f > min_f);

  std::mt19937_64 engine(seed);
  std::uniform_real_distribution<Real> random_dist(min_f, max_f);

  genRandArray(m_fm,      engine, random_dist);
  genRandArray(m_fvtheta, engine, random_dist);
  genRandArray(m_ft,      engine, random_dist);
  genRandArray(m_fphi,    engine, random_dist);
}

} // namespace Homme
