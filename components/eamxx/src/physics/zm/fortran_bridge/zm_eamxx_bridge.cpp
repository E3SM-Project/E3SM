#include "zm_eamxx_bridge.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

#include "physics/zm/fortran_bridge/zm_eamxx_bridge.hpp"

using scream::Real;
using scream::Int;

// A C++ interface to ZM fortran calls and vice versa

extern "C" {
  void zm_eamxx_bridge_init_c();
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init(){
  zm_eamxx_bridge_init_c();
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
