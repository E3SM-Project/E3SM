/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_REMAP_STATE_PROVIDER_HPP
#define HOMMEXX_REMAP_STATE_PROVIDER_HPP

#include "ElementsState.hpp"
#include "ErrorDefs.hpp"
#include "KernelVariables.hpp"
#include "Types.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme {
namespace Remap {

struct RemapStateProvider {

  ElementsState m_state;
  ExecViewManaged<Scalar* [NP][NP][NUM_LEV]> m_w_i_delta;
  ExecViewManaged<Scalar* [NP][NP][NUM_LEV]> m_phinh_i_delta;

  explicit RemapStateProvider(const ElementsState& state)
   : m_state(state)
   , m_w_i_delta ("w_i increments",state.num_elems())
   , m_phinh_i_delta ("phinh_i increments",state.num_elems())
  {}

  int num_states_remap() const { return 4;}

  void preprocess_states  (const int np1) const {
    // Compute w_i_delta and phinh_i_delta
  }

  void postprocess_states (const int np1) const {
    // Update w_i and phinh_i
  }

  KOKKOS_INLINE_FUNCTION
  ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>
  get_state(const KernelVariables &kv, int np1, int var) const {
    assert(var>=0 && var<=3);
    switch (var) {
    case 0:
      return Homme::subview(m_state.m_v, kv.ie, np1, 0);
    case 1:
      return Homme::subview(m_state.m_v, kv.ie, np1, 1);
    case 2:
      return Homme::subview(m_w_i_delta, kv.ie);
    case 3:
      return Homme::subview(m_phinh_i_delta, kv.ie);
    default:
      Errors::runtime_abort("RemapStateProvider: invalid variable index.\n");
      return ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>();
    }
  }
};

} // namespace Remap
} // namespace Homme

#endif // HOMMEXX_REMAP_STATE_PROVIDER_HPP
