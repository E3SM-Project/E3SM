/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "ElementsBuffers.hpp"

namespace Homme {

void ElementsBuffers::init(const int num_elems) {
  pressure    = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("Non-hydrostatic pressure", num_elems);
  pi          = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("Hydrostatic Pressure at midpoints", num_elems);
  pi_grad     = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]  >("Hydrostatic Pressure gradient at midpoints", num_elems);
  exner       = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("Exner pressure", num_elems);
  v_i         = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV_P]>("Horizontal velocity at interfaces", num_elems);
  dp_i        = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV_P]>("Pressure increments at interfaces", num_elems);
  dpnh_dp_i   = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV_P]>("Non-hydrostatic over hydrostatic pressure ratio", num_elems);
  vdp         = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]  >("(u,v)*dp", num_elems);
  div_vdp     = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("Divergence of (u,v)*dp", num_elems);
  ephi        = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("Energy (phi+0.5*|u|^2", num_elems);
  temp_column = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("A temporary column worth of data", num_elems);
  omega_p     = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("Omega_P = omega/pressure = (Dp/Dt)/pressure", num_elems);

  vstar       = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]  >("Buffer for (flux v)/dp", num_elems);
  dpdissk     = ExecViewManaged<Scalar *    [NP][NP][NUM_LEV]  >("dpdissk", num_elems);
  // pressure_grad = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Gradient of pressure", num_elems);
  // temperature_virt = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Virtual Temperature", num_elems);
  // temperature_grad = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Gradient of temperature", num_elems);
  // ephi        = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Kinetic Energy + Geopotential Energy", num_elems);
  // energy_grad = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Gradient of ephi", num_elems);
  // vorticity   = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Vorticity", num_elems);

  // ttens  = ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]>("Temporary for temperature",num_elems);
  // dptens = ExecViewManaged<Scalar*    [NP][NP][NUM_LEV]>("Temporary for dp3d",num_elems);
  // vtens  = ExecViewManaged<Scalar* [2][NP][NP][NUM_LEV]>("Temporary for velocity",num_elems);

  // vstar = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("buffer for (flux v)/dp",
  //      num_elems);
  // dpdissk = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>(
  //     "dpdissk", num_elems);

  // preq_buf = ExecViewManaged<Real * [NP][NP]>("Preq Buffer", num_elems);

  // sdot_sum = ExecViewManaged<Real * [NP][NP]>("Sdot sum buffer", num_elems);

  // div_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Divergence Buffer",
  //                                                          num_elems);
  // grad_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Gradient Buffer",
  //                                                           num_elems);
  // curl_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("Vorticity Buffer",
  //                                                           num_elems);

  // sphere_vector_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("laplacian vector Buffer", num_elems);

  // divergence_temp = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Divergence temporary",
  //                                                           num_elems);
  // vorticity_temp = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Vorticity temporary",
  //                                                           num_elems);
  // lapl_buf_1 = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Scalar laplacian Buffer", num_elems);
  // lapl_buf_2 = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Scalar laplacian Buffer", num_elems);
  // lapl_buf_3 = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("Scalar laplacian Buffer", num_elems);
  // v_vadv_buf = ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>("v_vadv buffer",
  //                                                             num_elems);
  // t_vadv_buf = ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>("t_vadv buffer",
  //                                                          num_elems);
  // eta_dot_dpdn_buf = ExecViewManaged<Scalar * [NP][NP][NUM_LEV_P]>("eta_dot_dpdpn buffer",
  //                                                                  num_elems);

  // kernel_start_times = ExecViewManaged<clock_t *>("Start Times", num_elems);
  // kernel_end_times = ExecViewManaged<clock_t *>("End Times", num_elems);
}

} // namespace Homme
