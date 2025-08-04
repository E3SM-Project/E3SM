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

  void zm_eamxx_bridge_run_c( Int  ncol,
                              Real dtime,
                              bool is_first_step,
                        const Real *state_phis,
                              Real *state_z_mid,
                              Real *state_z_int,
                        const Real *state_p_mid,
                        const Real *state_p_int,
                        const Real *state_p_del,
                              Real *state_t,
                              Real *state_qv,
                              Real *state_qc,
                              Real *state_u,
                              Real *state_v,
                              Real *state_omega,
                        const Real *state_cldfrac,
                        const Real *state_pblh,
                        const Real *landfrac,
                              Real *output_prec,
                              Real *output_snow,
                              Real *output_cape,
                              Real *output_tend_s,
                              Real *output_tend_q,
                              Real *output_tend_u,
                              Real *output_tend_v,
                              Real *output_rain_prod,
                              Real *output_snow_prod,
                              Real *output_prec_flux,
                              Real *output_snow_flux,
                              Real *output_mass_flux
                            );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcol, Int pver ){
  zm_eamxx_bridge_init_c( pcol, pver );
}

void zm_eamxx_bridge_run( Int ncol, Int pcol, Int pver,
                          ZMF::zm_input_state& zm_input,
                          ZMF::zm_output_tend& zm_output
                        ){
  //----------------------------------------------------------------------------
  zm_input.transpose<ekat::TransposeDirection::c2f>(ncol,pver);
  zm_output.transpose<ekat::TransposeDirection::c2f>(ncol,pver);

  // for (int i=0; i<zm_output.ncol; ++i) {
  //   for (int k=0; k<pver; ++k) {
  //     std::cout << "zm_eamxx_bridge_run (pre fortran) - ("<<i<<","<<k<<") - u / v : "
  //     << zm_input.f_u(i,k) << " / "
  //     << zm_input.f_v(i,k)
  //     << std::endl;
  //   }
  // }

  zm_eamxx_bridge_run_c( ncol,
                         zm_input.dtime,
                         zm_input.is_first_step,
                         zm_input.phis          .data(),
                         zm_input.f_z_mid       .data(),
                         zm_input.f_z_int       .data(),
                         zm_input.f_p_mid       .data(),
                         zm_input.f_p_int       .data(),
                         zm_input.f_p_del       .data(),
                         zm_input.f_T_mid       .data(),
                         zm_input.f_qv          .data(),
                         zm_input.f_qc          .data(),
                         zm_input.f_uwind       .data(),
                         zm_input.f_vwind       .data(),
                         zm_input.f_omega       .data(),
                         zm_input.f_cldfrac     .data(),
                         zm_input.pblh          .data(),
                         zm_input.landfrac      .data(),
                         zm_output.prec         .data(),
                         zm_output.snow         .data(),
                         zm_output.cape         .data(),
                         zm_output.f_tend_s     .data(),
                         zm_output.f_tend_qv    .data(),
                         zm_output.f_tend_u     .data(),
                         zm_output.f_tend_v     .data(),
                         zm_output.f_rain_prod  .data(),
                         zm_output.f_snow_prod  .data(),
                         zm_output.f_prec_flux  .data(),
                         zm_output.f_snow_flux  .data(),
                         zm_output.f_mass_flux  .data()
                        );

  zm_input.transpose<ekat::TransposeDirection::f2c>(ncol,pver);
  zm_output.transpose<ekat::TransposeDirection::f2c>(ncol,pver);
  //----------------------------------------------------------------------------
}

// end _c impls

} // namespace zm
} // namespace scream
