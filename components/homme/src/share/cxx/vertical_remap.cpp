/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"
#include "TimeLevel.hpp"

#include "VerticalRemapManager.hpp"

#include "Types.hpp"

namespace Homme
{

void vertical_remap(const Real dt)
{
  VerticalRemapManager& vrm = Context::singleton().get_vertical_remap_manager();
  TimeLevel& tl = Context::singleton().get_time_level();
  vrm.run_remap(tl.np1, tl.np1_qdp, dt);
}

} // namespace Homme
