/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/


#include <random>

#include "Tracers.hpp"
#include "Context.hpp"
#include "SimulationParams.hpp"

#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"

namespace Homme {

Tracers::Tracers(const int num_elems, const int num_tracers)
  : nt(num_tracers)
{
  Q = decltype(Q)("tracers concentration", num_elems);
  qdp = decltype(qdp)("tracers mass", num_elems);
  qtens_biharmonic = decltype(qtens_biharmonic)("qtens(_biharmonic)", num_elems);
  qlim = decltype(qlim)("qlim", num_elems);
}

void Tracers::random_init() {
  constexpr Real min_value = 0.015625;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(qdp, engine, random_dist);
  genRandArray(qtens_biharmonic, engine, random_dist);
  genRandArray(qlim, engine, random_dist);
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
