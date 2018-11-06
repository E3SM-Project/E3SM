/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_HPP
#define HOMMEXX_ELEMENTS_HPP

#include "Types.hpp"

namespace Homme {

class HybridVCoord;

/* Per element data - specific velocity, temperature, pressure, etc. */
class Elements {
public:
  // Coriolis term
  ExecViewManaged<Real * [NP][NP]> m_fcor;
  // Quadrature weights and metric tensor
  ExecViewManaged<Real * [NP][NP]>        m_mp;
  ExecViewManaged<Real * [NP][NP]>        m_spheremp;
  ExecViewManaged<Real * [NP][NP]>        m_rspheremp;
  ExecViewManaged<Real * [2][2][NP][NP]>  m_metinv;
  ExecViewManaged<Real * [NP][NP]>        m_metdet;
  ExecViewManaged<Real * [2][2][NP][NP]>  m_tensorvisc;
  ExecViewManaged<Real * [2][3][NP][NP]>  m_vec_sph2cart;
  // Prescribed surface geopotential height at eta = 1
  ExecViewManaged<Real * [NP][NP]> m_phis;

  // D (map for covariant coordinates) and D^{-1}
  ExecViewManaged<Real * [2][2][NP][NP]> m_d;
  ExecViewManaged<Real * [2][2][NP][NP]> m_dinv;

  // Omega is the 'pressure vertical velocity' in papers,
  // but omega=Dp/Dt  (not really vertical velocity).
  // In homme omega is scaled, derived%omega_p=(1/p)*(Dp/Dt)
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> m_omega_p;
  // Geopotential height field
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> m_phi;
  // weighted velocity flux for consistency
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]> m_derived_vn0;

  // Velocity in lon lat basis
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV]> m_v;
  // Temperature
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]> m_t;
  // dp ( it is dp/d\eta * delta(eta)), or pseudodensity
  ExecViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]> m_dp3d;

  ExecViewManaged<Real * [NUM_TIME_LEVELS][NP][NP]> m_ps_v;

  // eta=$\eta$ is the vertical coordinate
  // eta_dot_dpdn = $\dot{eta}\frac{dp}{d\eta}$
  //    (note there are NUM_PHYSICAL_LEV+1 of them)
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> m_eta_dot_dpdn;
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>
    m_derived_dp,                // for dp_tracers at physics timestep
    m_derived_divdp,             // divergence of dp
    m_derived_divdp_proj,        // DSSed divdp
    m_derived_dpdiss_biharmonic, // mean dp dissipation tendency, if nu_p>0
    m_derived_dpdiss_ave;        // mean dp used to compute psdiss_tens

  // Per Element Forcings
  // Momentum (? units are wrong in apply_cam_forcing...) forcing
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]> m_fm;
  // Temperature forcing
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> m_ft;

  //buffer views are temporaries that matter only during local RK steps
  //(dynamics and tracers time step).
  //m_ views are also used outside of local timesteps.
  struct BufferViews {

    BufferViews() = default;
    void init(const int num_elems);
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> pressure;
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> pressure_grad;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> temperature_virt;
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> temperature_grad;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> omega_p;
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> vdp;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> div_vdp;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> ephi;
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> energy_grad;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> vorticity;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> ttens;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> dptens;
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> vtens;

    // Buffers for EulerStepFunctor
    ExecViewManaged<Scalar*          [2][NP][NP][NUM_LEV]>  vstar;
    ExecViewManaged<Scalar*             [NP][NP][NUM_LEV]>  dpdissk;

    ExecViewManaged<Real* [NP][NP]> preq_buf;
    // sdot_sum is used in case rsplit=0 and in energy diagnostics
    // (not yet coded).
    ExecViewManaged<Real* [NP][NP]> sdot_sum;
    // Buffers for spherical operators
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> div_buf;
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> grad_buf;
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> curl_buf;

    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> sphere_vector_buf;

    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> divergence_temp;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> vorticity_temp;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> lapl_buf_1;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> lapl_buf_2;
    ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> lapl_buf_3;

    // Buffers for vertical advection terms in V and T for case
    // of Eulerian advection, rsplit=0. These buffers are used in both
    // cases, rsplit>0 and =0. Some of these values need to be init-ed
    // to zero at the beginning of each RK stage. Right now there is a code
    // for this, since Elements is a singleton.
    ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> v_vadv_buf;
    ExecViewManaged<Scalar* [NP][NP][NUM_LEV]> t_vadv_buf;
    ExecViewManaged<Scalar* [NP][NP][NUM_LEV_P]> eta_dot_dpdn_buf;

    ExecViewManaged<clock_t *> kernel_start_times;
    ExecViewManaged<clock_t *> kernel_end_times;
  } buffers;

  Elements() = default;

  void init(const int num_elems, const bool consthv);

  void random_init(int num_elems, Real max_pressure = 1.0);
  void random_init(int num_elems, Real max_pressure, const HybridVCoord& hvcoord);

  KOKKOS_INLINE_FUNCTION
  int num_elems() const { return m_num_elems; }

  // Fill the exec space views with data coming from F90 pointers
  void init_2d(const int ie, CF90Ptr &D, CF90Ptr &Dinv, CF90Ptr &fcor,
               CF90Ptr &mp, CF90Ptr &spheremp, CF90Ptr &rspheremp,
               CF90Ptr &metdet, CF90Ptr &metinv, 
               CF90Ptr &phis,
               CF90Ptr &tensorvisc,
               CF90Ptr &vec_sph2cart,
               const bool consthv);

  // Fill the exec space views with data coming from F90 pointers
  void pull_from_f90_pointers(CF90Ptr &state_v, CF90Ptr &state_t,
                              CF90Ptr &state_dp3d, CF90Ptr &derived_phi,
                              CF90Ptr &derived_omega_p,
                              CF90Ptr &derived_v, CF90Ptr &derived_eta_dot_dpdn);
  void pull_3d(CF90Ptr &derived_phi,
               CF90Ptr &derived_omega_p, CF90Ptr &derived_v);
  void pull_4d(CF90Ptr &state_v, CF90Ptr &state_t, CF90Ptr &state_dp3d);
  void pull_eta_dot(CF90Ptr &derived_eta_dot_dpdn);

  // Push the results from the exec space views to the F90 pointers
  void push_to_f90_pointers(F90Ptr &state_v, F90Ptr &state_t, F90Ptr &state_dp,
                            F90Ptr &derived_phi,
                            F90Ptr &derived_omega_p, F90Ptr &derived_v,
                            F90Ptr &derived_eta_dot_dpdn) const;
  void push_3d(F90Ptr &derived_phi,
               F90Ptr &derived_omega_p, F90Ptr &derived_v) const;
  void push_4d(F90Ptr &state_v, F90Ptr &state_t, F90Ptr &state_dp3d) const;
  void push_eta_dot(F90Ptr &derived_eta_dot_dpdn) const;

  void d(Real *d_ptr, int ie) const;
  void dinv(Real *dinv_ptr, int ie) const;

private:
  int m_num_elems;
};

} // Homme

#endif // HOMMEXX_ELEMENTS_HPP
