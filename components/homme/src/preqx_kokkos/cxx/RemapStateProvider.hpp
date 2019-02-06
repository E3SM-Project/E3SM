/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_REMAP_STATE_PROVIDER_HPP
#define HOMMEXX_REMAP_STATE_PROVIDER_HPP

#include "Types.hpp"
#include "ElementsState.hpp"
#include "ErrorDefs.hpp"
#include "KernelVariables.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme {
namespace Remap {

struct RemapStateProvider {

  ElementsState m_state;

  explicit RemapStateProvider(const ElementsState& state)
   : m_state(state) {}

  int num_states_remap() const { return 3;}

  void preprocess_states  (const int /* np1 */) const {
    // Nothing to do here for preqx
  }

  void postprocess_states (const int /* np1 */) const {
    // Nothing to do here for preqx
  }

  KOKKOS_INLINE_FUNCTION
  ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>
  get_state(const KernelVariables &kv, int np1, int var) const {
    assert(var>=0 && var<=2);
    switch (var) {
    case 0:
      return Homme::subview(m_state.m_v, kv.ie, np1, 0);
    case 1:
      return Homme::subview(m_state.m_v, kv.ie, np1, 1);
    case 2:
      return Homme::subview(m_state.m_t, kv.ie, np1);
    default:
      Errors::runtime_abort("RemapStateProvider: invalid variable index.\n");
      return ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>();
    }
  }
};

} // namespace Remap
} // namespace Homme

#endif // HOMMEXX_REMAP_STATE_PROVIDER_HPP
