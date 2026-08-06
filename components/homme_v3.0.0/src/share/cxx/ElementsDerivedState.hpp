/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_DERIVED_STATE_HPP
#define HOMMEXX_ELEMENTS_DERIVED_STATE_HPP

#include "Types.hpp"

namespace Homme {

// Per element derived data
class ElementsDerivedState {
public:

  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>     m_omega_p;  // Scaled 'pressure vertical velocity' (omega=(1/p)*Dp/Dt)
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>  m_vn0;      // weighted velocity flux for consistency
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>  m_vstar;    // velocity at start of tracer time step

  // eta=$\eta$ is the vertical coordinate
  // eta_dot_dpdn = $\dot{eta}\frac{dp}{d\eta}$
  // Note: the last level (surface) is always 0.
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV_P]>   m_eta_dot_dpdn;

  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>     m_dp;                // for dp_tracers at physics timestep
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>     m_divdp;             // divergence of dp
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>     m_divdp_proj;        // DSSed divdp
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>     m_dpdiss_biharmonic; // mean dp dissipation tendency, if nu_p>0
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>     m_dpdiss_ave;        // mean dp used to compute psdiss_tens

  ElementsDerivedState() : m_num_elems(0) {}

  void init (const int num_elems);

  void randomize(const int seed, const Real dp3d_min);

  KOKKOS_INLINE_FUNCTION
  int num_elems() const { return m_num_elems; }

private:
  int m_num_elems;
};

} // Homme

#endif // HOMMEXX_ELEMENTS_DERIVED_STATE_HPP
