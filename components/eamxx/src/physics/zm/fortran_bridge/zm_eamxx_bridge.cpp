#include "zm_eamxx_bridge.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

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
                              Real *state_qc,           // 12
                              Real *state_u,            // 13
                              Real *state_v,            // 14
                              Real *state_omega,        // 15
                        const Real *state_cldfrac,      // 16
                        const Real *state_pblh,         // 17
                        const Real *tpert,              // 18
                        const Real *landfrac,           // 19
                              Real *output_prec,        // 20
                              Real *output_snow,        // 21
                              Real *output_cape,        // 22
                              Int  *output_activity,    // 23
                              Real *output_tend_s,      // 24
                              Real *output_tend_q,      // 25
                              Real *output_tend_u,      // 26
                              Real *output_tend_v,      // 27
                              Real *output_rain_prod,   // 28
                              Real *output_snow_prod,   // 29
                              Real *output_prec_flux,   // 30
                              Real *output_snow_flux,   // 31
                              Real *output_mass_flux    // 32
                            );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcol, Int pver ){
  zm_eamxx_bridge_init_c( pcol, pver );
}

void zm_eamxx_bridge_run( Int ncol, Int pver,
                          ZMF::zm_input_state& zm_input,
                          ZMF::zm_output_tend& zm_output
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
                         zm_input.f_qc          .data(),  // 12
                         zm_input.f_uwind       .data(),  // 13
                         zm_input.f_vwind       .data(),  // 14
                         zm_input.f_omega       .data(),  // 15
                         zm_input.f_cldfrac     .data(),  // 16
                         zm_input.pblh          .data(),  // 17
                         zm_input.tpert         .data(),  // 18
                         zm_input.landfrac      .data(),  // 19
                         zm_output.prec         .data(),  // 20
                         zm_output.snow         .data(),  // 21
                         zm_output.cape         .data(),  // 22
                         zm_output.activity     .data(),  // 23
                         zm_output.f_tend_s     .data(),  // 24
                         zm_output.f_tend_qv    .data(),  // 25
                         zm_output.f_tend_u     .data(),  // 26
                         zm_output.f_tend_v     .data(),  // 27
                         zm_output.f_rain_prod  .data(),  // 28
                         zm_output.f_snow_prod  .data(),  // 29
                         zm_output.f_prec_flux  .data(),  // 30
                         zm_output.f_snow_flux  .data(),  // 31
                         zm_output.f_mass_flux  .data()   // 32
                        );

  zm_input.transpose<ekat::TransposeDirection::f2c>(ncol,pver);
  zm_output.transpose<ekat::TransposeDirection::f2c>(ncol,pver);

  //----------------------------------------------------------------------------
}

// end _c impls

} // namespace zm
} // namespace scream
