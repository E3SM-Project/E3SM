#include "zm_eamxx_bridge.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;

// A C++ interface to ZM fortran calls and vice versa

extern "C" {
  void zm_eamxx_bridge_init_c( Int pcols,
                               Int pver );
  void zm_eamxx_bridge_run_c( Int ncol,
                              Real *temperature,
                              Real *sp_humidity );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcols, Int pver ){
  std:cout << "zm_eamxx_bridge_init - pcols: " << pcols << std::endl;
  std:cout << "zm_eamxx_bridge_init - pver : " << pcols << std::endl;
  zm_eamxx_bridge_init_c( pver, pver );
}

void zm_eamxx_bridge_run(ZMF::zm_input_state& zm_input){
  zm_input.transpose<ekat::TransposeDirection::c2f>();
  zm_eamxx_bridge_run_c( zm_input.ncol,
                         ekat::scalarize(zm_input.T_mid).data(),
                         ekat::scalarize(zm_input.qv).data() );
}

// void zm_find_mse_max(zm_data_find_mse_max& d){
//   d.transpose<ekat::TransposeDirection::c2f>();
//   zm_find_mse_max_c( d.pcols,
//                      d.ncol,
//                      d.pver,
//                      d.num_msg,
//                      d.msemax_top_k,
//                      d.pergro_active,
//                      d.temperature,
//                      d.zmid,
//                      d.sp_humidity,
//                      d.msemax_klev,
//                      d.mse_max_val );
//   d.transpose<ekat::TransposeDirection::f2c>();
// }

// end _c impls

} // namespace zm
} // namespace scream
