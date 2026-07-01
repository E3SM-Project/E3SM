#include "zm_eamxx_bridge.hpp"

using scream::Real;
using scream::Int;

// A C++ interface to ZM fortran calls and vice versa

extern "C" {

void zm_eamxx_bridge_init_c( Int  pver_in, Int limcnv_in,
                             bool trig_dcape_in, bool trig_ull_in,
                             bool clos_dyn_adj_in, bool mcsp_enabled_in,
                             Real mcsp_t_coeff_in, Real mcsp_q_coeff_in,
                             Real mcsp_mom_coeff_in, bool mcsp_use_full_shear_in );

void zm_eamxx_bridge_run_c( Int  ncol,                // 01
                            Real dtime,               // 02
                            bool is_first_step,       // 03
                            const Real *state_phis,   // 04
                            Real *state_z_mid,        // 05
                            Real *state_z_int,        // 06
                            const Real *state_p_mid,  // 07
                            const Real *state_p_int,  // 08
                            const Real *state_p_del,  // 09
                            Real *state_t,            // 10
                            Real *state_qv,           // 11
                            Real *state_u,            // 12
                            Real *state_v,            // 13
                            Real *state_omega,        // 14
                            const Real *state_cldfrac,// 15
                            const Real *state_pblh,   // 16
                            const Real *tpert,        // 17
                            const Real *landfrac,     // 18
                            Real *t_star,             // 19 DCAPE T from time step n-1
                            Real *q_star,             // 20 DCAPE q from time step n-1
                            Real *output_prec,        // 21
                            Real *output_snow,        // 22
                            Real *output_cape,        // 23
                            Real *output_dcape,       // 24
                            Int  *output_activity,    // 25
                            Real *output_tend_t,      // 26
                            Real *output_tend_q,      // 27
                            Real *output_tend_u,      // 28
                            Real *output_tend_v,      // 29
                            Real *output_rain_prod,   // 30
                            Real *output_snow_prod,   // 31
                            Real *output_prec_flux,   // 32
                            Real *output_snow_flux,   // 33
                            Real *output_mass_flux,   // 34
                            Real *output_dlf,         // 35
                            Real *mcsp_freq,          // 36
                            Real *mcsp_shear,         // 37
                            Real *zm_depth,           // 38
                            Real *mcsp_ds_out,        // 39
                            Real *mcsp_dq_out,        // 40
                            Real *mcsp_du_out,        // 41
                            Real *mcsp_dv_out,        // 42
                            Real *evap_ds_out,        // 43
                            Real *evap_dq_out         // 44
                            );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pver, Int limcnv_in,
                           bool trig_dcape_in, bool trig_ull_in,
                           bool clos_dyn_adj_in, bool mcsp_enabled_in,
                           Real mcsp_t_coeff_in, Real mcsp_q_coeff_in,
                           Real mcsp_mom_coeff_in, bool mcsp_use_full_shear_in ){
  zm_eamxx_bridge_init_c( pver, limcnv_in,
                          trig_dcape_in, trig_ull_in,
                          clos_dyn_adj_in, mcsp_enabled_in,
                          mcsp_t_coeff_in, mcsp_q_coeff_in,
                          mcsp_mom_coeff_in, mcsp_use_full_shear_in );
}

void zm_eamxx_bridge_run( Int ncol, Int pver, Real dtime, bool is_first_step,
                          ZMF::ZmInputState& zm_input,
                          ZMF::ZmOutputTend& zm_output,
                          ZMF::ZmRuntimeOpt& zm_opts
                        ){
  //----------------------------------------------------------------------------
  zm_input.transpose<ekat::TransposeDirection::c2f>(ncol,pver);

  zm_eamxx_bridge_run_c( ncol,                            // 01
                         dtime,                           // 02
                         is_first_step,                   // 03
                         zm_input.h_phis        .data(),  // 04
                         zm_input.h_z_mid       .data(),  // 05
                         zm_input.h_z_int       .data(),  // 06
                         zm_input.h_p_mid       .data(),  // 07
                         zm_input.h_p_int       .data(),  // 08
                         zm_input.h_p_del       .data(),  // 09
                         zm_input.h_T_mid       .data(),  // 10
                         zm_input.h_qv          .data(),  // 11
                         zm_input.h_uwind       .data(),  // 12
                         zm_input.h_vwind       .data(),  // 13
                         zm_input.h_omega       .data(),  // 14
                         zm_input.h_cldfrac     .data(),  // 15
                         zm_input.h_pblh        .data(),  // 16
                         zm_input.h_tpert       .data(),  // 17
                         zm_input.h_landfrac    .data(),  // 18
                         zm_input.h_t_prev      .data(),  // 19
                         zm_input.h_q_prev      .data(),  // 20
                         zm_output.h_prec       .data(),  // 21
                         zm_output.h_snow       .data(),  // 22
                         zm_output.h_cape       .data(),  // 23
                         zm_output.h_dcape      .data(),  // 24
                         zm_output.h_activity   .data(),  // 25
                         zm_output.h_tend_t     .data(),  // 26
                         zm_output.h_tend_qv    .data(),  // 27
                         zm_output.h_tend_u     .data(),  // 28
                         zm_output.h_tend_v     .data(),  // 29
                         zm_output.h_rain_prod  .data(),  // 30
                         zm_output.h_snow_prod  .data(),  // 31
                         zm_output.h_prec_flux  .data(),  // 32
                         zm_output.h_snow_flux  .data(),  // 33
                         zm_output.h_mass_flux  .data(),  // 34
                         zm_output.h_dlf        .data(),  // 35
                         zm_output.h_mcsp_freq  .data(),  // 36
                         zm_output.h_mcsp_shear .data(),  // 37
                         zm_output.h_zm_depth   .data(),  // 38
                         zm_output.h_mcsp_ds_out.data(),  // 39
                         zm_output.h_mcsp_dq_out.data(),  // 40
                         zm_output.h_mcsp_du_out.data(),  // 41
                         zm_output.h_mcsp_dv_out.data(),  // 42
                         zm_output.h_evap_ds_out.data(),  // 43
                         zm_output.h_evap_dq_out.data()   // 44
                        );

  zm_output.transpose<ekat::TransposeDirection::f2c>(ncol,pver);

  //----------------------------------------------------------------------------
}

// end _c impls

} // namespace zm
} // namespace scream
