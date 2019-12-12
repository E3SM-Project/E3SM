/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_REMAP_STATE_PROVIDER_HPP
#define HOMMEXX_REMAP_STATE_PROVIDER_HPP

#include "Types.hpp"
#include "Elements.hpp"
#include "ErrorDefs.hpp"
#include "KernelVariables.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme {
namespace Remap {

struct RemapStateProvider {

  ElementsState m_state;

  explicit RemapStateProvider(const Elements& elements)
   : m_state(elements.m_state) {}

  template<typename Tag>
  void allocate_buffers(const Kokkos::TeamPolicy<ExecSpace,Tag>& /* policy */) {
    // preqx does not do anything during states pre/post processing,
    // so no need to create any buffer
  }

  KOKKOS_INLINE_FUNCTION
  int num_states_remap() const { return 3;}

  KOKKOS_INLINE_FUNCTION
  int num_states_preprocess() const { return 0;}

  KOKKOS_INLINE_FUNCTION
  int num_states_postprocess() const { return 0;}

  KOKKOS_INLINE_FUNCTION
  bool is_intrinsic_state (const int istate) const {
    // Sanity check
    assert (istate>=0 && istate<num_states_remap());
    static_assert(sizeof(decltype(istate))>0,"Silence compiler warning went wrong.\n");

    // All states are intrinsic in preqx, so they all need to be multiplied by dp
    // during the remap procedure

    return true;
  }

  KOKKOS_INLINE_FUNCTION
  void preprocess_state (const KernelVariables& /* kv */,
                         const int /* istate */,
                         const int /* np1 */,
                         ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> /* dp */) const {
    // Nothing to do here for preqx
  }

  KOKKOS_INLINE_FUNCTION
  void postprocess_state (const KernelVariables& /* kv */,
                          const int /* istate */,
                          const int /* np1 */,
                          ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> /* dp */) const {
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
      Kokkos::abort("RemapStateProvider: invalid variable index.\n");
      return ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>();
    }
  }
};

} // namespace Remap
} // namespace Homme

#endif // HOMMEXX_REMAP_STATE_PROVIDER_HPP
