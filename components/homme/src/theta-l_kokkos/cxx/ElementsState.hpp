/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_STATE_HPP
#define HOMMEXX_ELEMENTS_STATE_HPP

#include "Types.hpp"
#include "kokkos_utils.hpp"

namespace Homme {

class HybridVCoord;

// Reference states, needed in HV and vert remap
struct RefStates {
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV_P]> phi_i_ref;
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV  ]> theta_ref;
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV  ]> dp_ref;

  RefStates () :
    m_num_elems(0)
    , m_policy(get_default_team_policy<ExecSpace>(1))
    , m_tu(m_policy)
  {}

  void init (const int num_elems);

  int num_elems () const { return m_num_elems; }
private:
  int m_num_elems;
  Kokkos::TeamPolicy<ExecSpace> m_policy;
  TeamUtils<ExecSpace> m_tu;
};

/* Per element data - specific velocity, temperature, pressure, etc. */
class ElementsState {
public:

  RefStates m_ref_states;

  ExecViewManaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV  ]> m_v;          // Horizontal velocity
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV_P]> m_w_i;        // Vertical velocity at interfaces
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV  ]> m_vtheta_dp;  // Virtual potential temperature (mass)
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV_P]> m_phinh_i;    // Geopotential used by NH model at interfaces
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS]   [NP][NP][NUM_LEV  ]> m_dp3d;       // Delta p on levels

  ExecViewManaged<Real   * [NUM_TIME_LEVELS]   [NP][NP]           > m_ps_v;       // Surface pressure

  void init_storage(const int num_elems);

  ElementsState() :
    m_num_elems(0)
    , m_policy(get_default_team_policy<ExecSpace>(1))
    , m_tu(m_policy)
  {}

  void init(const int num_elems);

  void randomize(const int seed);
  void randomize(const int seed, const Real max_pressure);
  void randomize(const int seed, const Real max_pressure, const Real ps0, const Real hyai0);

  void randomize(const int seed, const Real max_pressure, const Real ps0, const Real hyai0,
                 const ExecViewUnmanaged<const Real*[NP][NP]>& phis);

  KOKKOS_INLINE_FUNCTION
  int num_elems() const { return m_num_elems; }

  // Fill the exec space views with data coming from F90 pointers
  void pull_from_f90_pointers(CF90Ptr& state_v,         CF90Ptr& state_w_i,
                              CF90Ptr& state_vtheta_dp, CF90Ptr& state_phinh_i,
                              CF90Ptr& state_dp3d,      CF90Ptr& state_ps_v);

  // Push the results from the exec space views to the F90 pointers
  void push_to_f90_pointers(F90Ptr& state_v, F90Ptr& state_w_i, F90Ptr& state_vtheta_dp,
                            F90Ptr& state_phinh_i, F90Ptr& state_dp) const;

private:
  int m_num_elems;
  Kokkos::TeamPolicy<ExecSpace> m_policy;
  TeamUtils<ExecSpace> m_tu;
};

// Check ElementsState for NaN or incorrectly signed values. The initial check
// is fast and on device. If everything is fine, the routine returns
// immediately. If there is a bad value, a subsequent check is run on the host,
// and this check prints detailed information to a file called
// hommexx.errlog.${rank}. Then runtime_abort is called with a message pointing
// to this file.
void check_print_abort_on_bad_elems(const std::string& label,    // string to ID call site
                                    const int time_level,        // time level index in state arrays
                                    const int error_code = -1);

} // Homme

#endif // HOMMEXX_ELEMENTS_STATE_HPP
