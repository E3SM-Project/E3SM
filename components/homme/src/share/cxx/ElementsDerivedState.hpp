/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_DERIVED_STATE_HPP
#define HOMMEXX_ELEMENTS_DERIVED_STATE_HPP

#include "Types.hpp"

namespace Homme {

/* Per element data - specific velocity, temperature, pressure, etc. */
class ElementsDerivedState {
public:

  // Omega is the 'pressure vertical velocity' in papers,
  // but omega=Dp/Dt  (not really vertical velocity).
  // In homme omega is scaled, derived%omega_p=(1/p)*(Dp/Dt)
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> m_omega_p;
  // Geopotential height field
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> m_phi;
  // weighted velocity flux for consistency
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]> m_vn0;

  // eta=$\eta$ is the vertical coordinate
  // eta_dot_dpdn = $\dot{eta}\frac{dp}{d\eta}$
  //    (note there are NUM_PHYSICAL_LEV+1 of them)
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>   m_eta_dot_dpdn;

  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>   m_dp;                // for dp_tracers at physics timestep
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>   m_divdp;             // divergence of dp
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>   m_divdp_proj;        // DSSed divdp
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>   m_dpdiss_biharmonic; // mean dp dissipation tendency, if nu_p>0
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>   m_dpdiss_ave;        // mean dp used to compute psdiss_tens

  // Per Element Forcings
  // Momentum (? units are wrong in apply_cam_forcing...) forcing
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]> m_fm;
  // Temperature forcing
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> m_ft;

  ElementsDerivedState() = default;

  void init(const int num_elems);

  void random_init(int num_elems, Real dp3d_min);

  KOKKOS_INLINE_FUNCTION
  int num_elems() const { return m_num_elems; }

  // Fill the exec space views with data coming from F90 pointers
  void pull_from_f90_pointers(CF90Ptr& phi, CF90Ptr& omega_p,
                              CF90Ptr& vn0, CF90Ptr& eta_dot_dpdn);
  void pull_3d(CF90Ptr& phi, CF90Ptr& omega_p, CF90Ptr& vn0);
  void pull_eta_dot(CF90Ptr &derived_eta_dot_dpdn);

  // Push the results from the exec space views to the F90 pointers
  void push_to_f90_pointers(F90Ptr& phi, F90Ptr& omega_p,
                            F90Ptr &vn0, F90Ptr& eta_dot_dpdn) const;
  void push_3d(F90Ptr& phi, F90Ptr& omega_p, F90Ptr& vn0) const;
  void push_eta_dot(F90Ptr &derived_eta_dot_dpdn) const;

private:
  int m_num_elems;
};

} // Homme

#endif // HOMMEXX_ELEMENTS_DERIVED_STATE_HPP
