#include "zm_eamxx_bridge.hpp"

using scream::Real;
using scream::Int;

// A C++ interface to ZM fortran calls and vice versa

extern "C" {
  void zm_eamxx_bridge_init_c(Int  pcol_in,
                              Int  pver_in );

  void zm_eamxx_bridge_run_c( Int  ncol,                // 01
                              Real dtime,               // 02
                              bool is_first_step,       // 03
                        const Real *state_phis,         // 04
                              Real *state_z_mid,        // 05
                              Real *state_z_int,        // 06
                        const Real *state_p_mid,        // 07
                        const Real *state_p_int,        // 08
                        const Real *state_p_del,        // 09
                              Real *state_t,            // 10
                              Real *state_qv,           // 11
                              Real *state_u,            // 12
                              Real *state_v,            // 13
                              Real *state_omega,        // 14
                        const Real *state_cldfrac,      // 15
                        const Real *state_pblh,         // 16
                        const Real *tpert,              // 17
                        const Real *landfrac,           // 18
                              Real *output_prec,        // 19
                              Real *output_snow,        // 20
                              Real *output_cape,        // 21
                              Int  *output_activity,    // 22
                              Real *output_tend_s,      // 23
                              Real *output_tend_q,      // 24
                              Real *output_tend_u,      // 25
                              Real *output_tend_v,      // 26
                              Real *output_rain_prod,   // 27
                              Real *output_snow_prod,   // 28
                              Real *output_prec_flux,   // 29
                              Real *output_snow_flux,   // 30
                              Real *output_mass_flux    // 31
                            );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcol, Int pver ){
  zm_eamxx_bridge_init_c( pcol, pver );
}

void zm_eamxx_bridge_run( Int ncol, Int pver,
                          ZMF::zm_input_state& zm_input,
                          ZMF::zm_output_tend& zm_output,
                          ZMF::zm_runtime_opt& zm_opts
                        ){
  //----------------------------------------------------------------------------
  zm_input.transpose<ekat::TransposeDirection::c2f>(ncol,pver);
  zm_output.transpose<ekat::TransposeDirection::c2f>(ncol,pver);

  zm_eamxx_bridge_run_c( ncol,                            // 01
                         zm_input.dtime,                  // 02
                         zm_input.is_first_step,          // 03
                         zm_input.phis          .data(),  // 04
                         zm_input.f_z_mid       .data(),  // 05
                         zm_input.f_z_int       .data(),  // 06
                         zm_input.f_p_mid       .data(),  // 07
                         zm_input.f_p_int       .data(),  // 08
                         zm_input.f_p_del       .data(),  // 09
                         zm_input.f_T_mid       .data(),  // 10
                         zm_input.f_qv          .data(),  // 11
                         zm_input.f_uwind       .data(),  // 12
                         zm_input.f_vwind       .data(),  // 13
                         zm_input.f_omega       .data(),  // 14
                         zm_input.f_cldfrac     .data(),  // 15
                         zm_input.pblh          .data(),  // 16
                         zm_input.tpert         .data(),  // 17
                         zm_input.landfrac      .data(),  // 18
                         zm_output.prec         .data(),  // 19
                         zm_output.snow         .data(),  // 20
                         zm_output.cape         .data(),  // 21
                         zm_output.activity     .data(),  // 22
                         zm_output.f_tend_s     .data(),  // 23
                         zm_output.f_tend_qv    .data(),  // 24
                         zm_output.f_tend_u     .data(),  // 25
                         zm_output.f_tend_v     .data(),  // 26
                         zm_output.f_rain_prod  .data(),  // 27
                         zm_output.f_snow_prod  .data(),  // 28
                         zm_output.f_prec_flux  .data(),  // 29
                         zm_output.f_snow_flux  .data(),  // 30
                         zm_output.f_mass_flux  .data()   // 31
                        );

  zm_input.transpose<ekat::TransposeDirection::f2c>(ncol,pver);
  zm_output.transpose<ekat::TransposeDirection::f2c>(ncol,pver);

  //----------------------------------------------------------------------------
}

// end _c impls

} // namespace zm
} // namespace scream
