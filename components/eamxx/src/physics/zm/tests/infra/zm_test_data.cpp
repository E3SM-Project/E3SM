#include "zm_test_data.hpp"

#include <ekat_pack_kokkos.hpp>

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to ZM fortran calls and vice versa
//

namespace scream {
namespace zm {

extern "C" {

  void zm_find_mse_max_c( Int  pcols,
                          Int  ncol,
                          Int  pver,
                          Int  num_msg,
                          Int  *msemax_top_k,
                          bool pergro_active,
                          Real *temperature,
                          Real *zmid,
                          Real *sp_humidity,
                          Int *msemax_klev,
                          Real *mse_max_val );

} // extern "C" : end _c decls

void zm_find_mse_max(zm_data_find_mse_max& d){
  d.transition<ekat::TransposeDirection::c2f>();
  zm_find_mse_max_c( d.pcols,
                     d.ncol,
                     d.pver,
                     d.num_msg,
                     d.msemax_top_k,
                     d.pergro_active,
                     d.temperature,
                     d.zmid,
                     d.sp_humidity,
                     d.msemax_klev,
                     d.mse_max_val );
  d.transition<ekat::TransposeDirection::f2c>();
}

// end _c impls

} // namespace zm
} // namespace scream
