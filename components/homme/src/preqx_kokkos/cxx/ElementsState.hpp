/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_STATE_HPP
#define HOMMEXX_ELEMENTS_STATE_HPP

#include "Types.hpp"
#include "utilities/Hash.hpp"

namespace Homme {

class HybridVCoord;

/* Per element data - specific velocity, temperature, pressure, etc. */
class ElementsState {
public:
  // Velocity in lon lat basis
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV]> m_v;
  // Temperature
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]> m_t;
  // dp ( it is dp/d\eta * delta(eta)), or pseudodensity
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]> m_dp3d;

  ExecViewManaged<Real * [NUM_TIME_LEVELS][NP][NP]> m_ps_v;

  ElementsState() : m_num_elems(0) {}

  void init(const int num_elems);

  void randomize(const int seed);
  void randomize(const int seed, const Real max_pressure);
  void randomize(const int seed, const Real max_pressure, const Real ps0, const Real hyai0);

  // phis is not used in preqx model initialization
  void randomize(const int seed, const Real max_pressure, const Real ps0, const Real hyai0,
                 const ExecViewUnmanaged<const Real*[NP][NP]>& /* phis */) {
    randomize(seed, max_pressure, ps0, hyai0);
  }

  KOKKOS_INLINE_FUNCTION
  int num_elems() const { return m_num_elems; }

  // Fill the exec space views with data coming from F90 pointers
  void pull_from_f90_pointers(CF90Ptr& state_v,    CF90Ptr& state_t,
                              CF90Ptr& state_dp3d, CF90Ptr& state_ps_v);

  // Push the results from the exec space views to the F90 pointers
  void push_to_f90_pointers(F90Ptr& state_v, F90Ptr& state_t, F90Ptr& state_dp) const;

  HashType hash(const int time_level) const;

private:
  int m_num_elems;
};

// Not implemented.
void check_print_abort_on_bad_elems(const std::string& label, const int time_level,
                                    const int error_code = -1);

} // Homme

#endif // HOMMEXX_ELEMENTS_STATE_HPP
