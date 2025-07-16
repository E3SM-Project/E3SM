#include "zm_eamxx_bridge.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;

// A C++ interface to ZM fortran calls and vice versa

extern "C" {
  void zm_eamxx_bridge_init_c(Int  pcols_in,
                              Int  pver_in );

  void zm_eamxx_bridge_run_c( Int  ncol,
                              bool is_first_step,
                              // Real *state_phis,
                              const Real *state_p_mid,
                              const Real *state_p_int,
                              const Real *state_p_del,
                              Real *state_t,
                              Real *state_qv,
                              Real *state_qc );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcols, Int pver ){
  zm_eamxx_bridge_init_c( pcols, pver );
}

void zm_eamxx_bridge_run(ZMF::zm_input_state& zm_input){

  // int ncol = zm_input.T_mid.extent(0);
  // int nlev = zm_input.T_mid.extent(1);

  //----------------------------------------------------------------------------
  auto s_p = ekat::scalarize( zm_input.p_mid );
  auto s_T = ekat::scalarize( zm_input.T_mid );
  auto s_q = ekat::scalarize( zm_input.qv    );
  for (int i=0; i<zm_input.ncol; ++i) {
    for (int k=0; k<12; ++k) {
      std::cout << "zm_eamxx_bridge_run - p("<<i<<","<<k<<"): " << s_p(i,k) << std::endl;
    }
  }

  zm_input.transpose<ekat::TransposeDirection::c2f>();

  zm_eamxx_bridge_run_c( zm_input.ncol,
                         zm_input.is_first_step,
                         // zm_input.phis    .data(),
                         zm_input.f_p_mid .data(),
                         zm_input.f_p_int .data(),
                         zm_input.f_p_del .data(),
                         zm_input.f_T_mid .data(),
                         zm_input.f_qv    .data(),
                         zm_input.f_qc    .data() );
  zm_input.transpose<ekat::TransposeDirection::f2c>(pver);

  // zm_eamxx_bridge_run_c( zm_input.ncol,
  //                        zm_input.is_first_step,
  //                        zm_input.phis  .data(),
  //                        zm_input.T_mid .data(),
  //                        zm_input.qv    .data() );
  
  // zm_output.transpose<ekat::TransposeDirection::f2c>();
}

// end _c impls

} // namespace zm
} // namespace scream
