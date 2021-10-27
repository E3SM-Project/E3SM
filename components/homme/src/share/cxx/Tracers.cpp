/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/


#include <random>

#include "Tracers.hpp"
#include "Context.hpp"
#include "Elements.hpp"
#include "SimulationParams.hpp"

#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"

namespace Homme {

Tracers::Tracers(const int num_elems, const int num_tracers)
{
  init (num_elems, num_tracers);
}

void Tracers::init(const int num_elems, const int num_tracers)
{
  // Sanity check
  assert(num_elems>0);
  assert(num_tracers>=0);

  ne = num_elems;
  nt = num_tracers;

  qdp = decltype(qdp)("tracers mass", num_elems);
  qtens_biharmonic = decltype(qtens_biharmonic)("qtens(_biharmonic)", num_elems);
  qlim = decltype(qlim)("qlim", num_elems);

  Q = decltype(Q)("tracers concentration", num_elems,num_tracers);
  fq = decltype(fq)("fq",num_elems,num_tracers);

  m_inited = true;
}

void Tracers::randomize(const int seed) {
  // Check tracers were inited
  assert (m_inited);

  constexpr Real min_value = 0.015625;
  std::mt19937_64 engine(seed);
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(qdp, engine, random_dist);
  genRandArray(qtens_biharmonic, engine, random_dist);
  genRandArray(qlim, engine, random_dist);
  genRandArray(fq, engine, random_dist);
  genRandArray(Q, engine, random_dist);
}

void Tracers::pull_qdp(CF90Ptr &state_qdp) {
  HostViewUnmanaged<
      const Real * [Q_NUM_TIME_LEVELS][QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]>
  state_qdp_f90(state_qdp, qdp.extent_int(0));
  sync_to_device(state_qdp_f90, qdp);
}

void Tracers::push_qdp(F90Ptr &state_qdp) const {
  HostViewUnmanaged<
      Real * [Q_NUM_TIME_LEVELS][QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]>
  state_qdp_f90(state_qdp, qdp.extent_int(0));
  sync_to_host(qdp, state_qdp_f90);
}

} // namespace Homme
