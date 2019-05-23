/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"
#include "CamForcing.hpp"
#include "Tracers.hpp"
#include "Elements.hpp"
#include "TimeLevel.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "KernelVariables.hpp"
#include "vector/vector_pragmas.hpp"
#include "profiling.hpp"

namespace Homme {


void apply_cam_forcing(const Real &dt) {
  GPTLstart("ApplyCAMForcing");
  GPTLstop("ApplyCAMForcing");
}

void apply_cam_forcing_dynamics(const Real &dt) {
  GPTLstart("ApplyCAMForcing_dynamics");
  GPTLstop("ApplyCAMForcing_dynamics");
}

} // namespace Homme
