/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"
#include "CamForcing.hpp"
#include "ForcingFunctor.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "profiling.hpp"

#include "utilities/MathUtils.hpp"

namespace Homme {

void apply_cam_forcing(const Real dt) {
  GPTLstart("ApplyCAMForcing");

  const auto& p  = Context::singleton().get<SimulationParams>();
  const auto& tl = Context::singleton().get<TimeLevel>();
  auto& ff = Context::singleton().get<ForcingFunctor>();
 
  ff.tracers_forcing(dt,tl.n0,tl.n0_qdp,false,p.moisture);
  ff.states_forcing(dt,tl.n0);

  GPTLstop("ApplyCAMForcing");
}

void apply_cam_forcing_dynamics(const Real dt) {
  GPTLstart("ApplyCAMForcing_dynamics");

  const auto& tl = Context::singleton().get<TimeLevel>();
  auto& ff = Context::singleton().get<ForcingFunctor>();

  ff.states_forcing(dt,tl.np1);

  GPTLstop("ApplyCAMForcing_dynamics");
}

} // namespace Homme
