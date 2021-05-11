/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ELEMENTS_BUFFERS_HPP
#define HOMMEXX_ELEMENTS_BUFFERS_HPP

#include "ElementsBuffersBase.hpp"


namespace Homme {

// Elements-dependent buffers.
struct ElementsBuffers : public ElementsBuffersBase {

  ElementsBuffers() = default;
  void init(const int num_elems);

  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   pressure;
  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   pi;
  ExecViewUnmanaged<Scalar* [2][NP][NP][NUM_LEV]  >   pi_grad;
  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   exner;
  ExecViewUnmanaged<Scalar* [2][NP][NP][NUM_LEV_P]>   v_i;
  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV_P]>   dp_i;
  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV_P]>   dpnh_dp_i;
  ExecViewUnmanaged<Scalar* [2][NP][NP][NUM_LEV]  >   vdp;
  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   div_vdp;
  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   ephi;
  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   omega_p;

  ExecViewUnmanaged<Scalar*    [NP][NP][NUM_LEV]  >   temp_column;

  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> pressure_grad;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> temperature_virt;
  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> temperature_grad;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> omega_p;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> ephi;
  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> energy_grad;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> vorticity;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> ttens;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> dptens;
  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> vtens;


  // Buffers for EulerStepFunctor
  ExecViewUnmanaged<Scalar*          [2][NP][NP][NUM_LEV]>  vstar;
  ExecViewUnmanaged<Scalar*             [NP][NP][NUM_LEV]>  dpdissk;

  // ExecViewManaged<Real* [NP][NP]> preq_buf;
  // // sdot_sum is used in case rsplit=0 and in energy diagnostics
  // // (not yet coded).
  // ExecViewManaged<Real* [NP][NP]> sdot_sum;
  // // Buffers for spherical operators
  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> div_buf;
  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> grad_buf;
  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> curl_buf;

  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> sphere_vector_buf;

  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> divergence_temp;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> vorticity_temp;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> lapl_buf_1;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> lapl_buf_2;
  // ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]> lapl_buf_3;

  // // Buffers for vertical advection terms in V and T for case
  // // of Eulerian advection, rsplit=0. These buffers are used in both
  // // cases, rsplit>0 and =0. Some of these values need to be init-ed
  // // to zero at the beginning of each RK stage.
  // ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]> v_vadv_buf;
  // ExecViewManaged<Scalar* [NP][NP][NUM_LEV]> t_vadv_buf;
  // ExecViewManaged<Scalar* [NP][NP][NUM_LEV_P]> eta_dot_dpdn_buf;

  // ExecViewManaged<clock_t *> kernel_start_times;
  // ExecViewManaged<clock_t *> kernel_end_times;
};

} // Homme

#endif // HOMMEXX_ELEMENTS_BUFFERS_HPP
