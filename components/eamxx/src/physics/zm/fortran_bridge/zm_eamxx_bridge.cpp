#include "zm_eamxx_bridge.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;

// A C++ interface to ZM fortran calls and vice versa

extern "C" {
  void zm_eamxx_bridge_init_c(Int  &pcols_in,
                              Int  &pver_in );

  void zm_eamxx_bridge_run_c( Int  &ncol,
                              bool &is_first_step,
                              Real *state_t,
                              Real *state_q );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcols, Int pver ){
  zm_eamxx_bridge_init_c( pcols, pver );
}

void zm_eamxx_bridge_run(ZMF::zm_input_state& zm_input){
  zm_input.transpose<ekat::TransposeDirection::c2f>();
  zm_eamxx_bridge_run_c( zm_input.ncol,
                         zm_input.is_first_step,
                         ekat::scalarize(zm_input.T_mid ).data(),
                         ekat::scalarize(zm_input.qv    ).data() );
  // zm_output.transpose<ekat::TransposeDirection::f2c>();
}

// end _c impls

} // namespace zm
} // namespace scream
