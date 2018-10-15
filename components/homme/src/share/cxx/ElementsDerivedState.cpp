/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "ElementsDerivedState.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"

#include <limits>
#include <random>

namespace Homme {

void ElementsDerivedState::init(const int num_elems) {
  m_num_elems = num_elems;

  m_omega_p = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Omega P", m_num_elems);
  m_phi     = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("PHI", m_num_elems);

  m_vn0 = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Derived Lateral Velocities", m_num_elems);

  m_fm = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("F_Momentum", m_num_elems);
  m_ft = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]>("F_Temperature", m_num_elems);

  m_eta_dot_dpdn = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("eta_dot_dpdn", m_num_elems);

  m_dp                = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dp", m_num_elems);
  m_divdp             = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_divdp", m_num_elems);
  m_divdp_proj        = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_divdp_proj", m_num_elems);
  m_dpdiss_biharmonic = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dpdiss_biharmonic", m_num_elems);
  m_dpdiss_ave        = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("derived_dpdiss_ave", m_num_elems);
}

//test for tensor hv is needed
void ElementsDerivedState::random_init(int num_elems, Real dp3d_min) {
  init(num_elems);

  // arbitrary minimum value to generate and minimum determinant allowed
  constexpr const Real min_value = 0.015625;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<Real> random_dist(min_value, 1.0 / min_value);

  genRandArray(m_omega_p, engine, random_dist);
  genRandArray(m_phi,     engine, random_dist);
  genRandArray(m_vn0,     engine, random_dist);

  // Generate eta_dot_dpdn so that it is << dp3d
  genRandArray(m_eta_dot_dpdn, engine, std::uniform_real_distribution<Real>(0.01*dp3d_min,0.1*dp3d_min));
}

void ElementsDerivedState::pull_from_f90_pointers(
    CF90Ptr& phi, CF90Ptr& omega_p,
    CF90Ptr& vn0, CF90Ptr& eta_dot_dpdn) {
  pull_3d(phi, omega_p, vn0);
  pull_eta_dot(eta_dot_dpdn);
}

void ElementsDerivedState::pull_3d(CF90Ptr& phi, CF90Ptr& omega_p, CF90Ptr& vn0) {
  HostViewUnmanaged<const Real *[NUM_PHYSICAL_LEV]   [NP][NP]> phi_f90(phi,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_PHYSICAL_LEV]   [NP][NP]> omega_p_f90(omega_p,m_num_elems);
  HostViewUnmanaged<const Real *[NUM_PHYSICAL_LEV][2][NP][NP]> vn0_f90(vn0,m_num_elems);

  sync_to_device(phi_f90,     m_phi);
  sync_to_device(omega_p_f90, m_omega_p);
  sync_to_device(vn0_f90,       m_vn0);
}

void ElementsDerivedState::pull_eta_dot(CF90Ptr& eta_dot_dpdn) {
  HostViewUnmanaged<const Real *[NUM_INTERFACE_LEV][NP][NP]> eta_dot_dpdn_f90(eta_dot_dpdn,m_num_elems);
  sync_to_device_i2p(eta_dot_dpdn_f90,m_eta_dot_dpdn);
}

void ElementsDerivedState::push_to_f90_pointers(F90Ptr& phi, F90Ptr& omega_p,
                                                F90Ptr& vn0, F90Ptr& eta_dot_dpdn) const {
  push_3d(phi, omega_p, vn0);
  push_eta_dot(eta_dot_dpdn);
}

void ElementsDerivedState::push_3d(F90Ptr& phi, F90Ptr& omega_p, F90Ptr& vn0) const {
  HostViewUnmanaged<Real *[NUM_PHYSICAL_LEV]   [NP][NP]> phi_f90    (phi,m_num_elems);
  HostViewUnmanaged<Real *[NUM_PHYSICAL_LEV]   [NP][NP]> omega_p_f90(omega_p,m_num_elems);
  HostViewUnmanaged<Real *[NUM_PHYSICAL_LEV][2][NP][NP]> vn0_f90    (vn0,m_num_elems);

  sync_to_host(m_phi,     phi_f90);
  sync_to_host(m_omega_p, omega_p_f90);
  sync_to_host(m_vn0,     vn0_f90);
}

void ElementsDerivedState::push_eta_dot(F90Ptr& eta_dot_dpdn) const {
  HostViewUnmanaged<Real *[NUM_INTERFACE_LEV][NP][NP]> eta_dot_dpdn_f90(eta_dot_dpdn,m_num_elems);
  sync_to_host_p2i(m_eta_dot_dpdn,eta_dot_dpdn_f90);
}

} // namespace Homme
