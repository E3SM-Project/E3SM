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
                              // Real *state_t,
                              // Real *state_q );
                              Real *state_phis,
                              Real *state_t,
                              Real *state_q );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcols, Int pver ){
  zm_eamxx_bridge_init_c( pcols, pver );
}

void zm_eamxx_bridge_run(ZMF::zm_input_state& zm_input){

  // int ncol = zm_input.T_mid.extent(0);
  // int nlev = zm_input.T_mid.extent(1);

  auto phis_s = zm_input.phis;
  auto T_s    = ekat::scalarize( zm_input.T_mid );
  auto q_s    = ekat::scalarize( zm_input.qv    );

  int ncol = T_s.extent(0);
  int nlev = T_s.extent(1);

  std::cout << "zm_eamxx_bridge_run - ncol: "<< ncol << std::endl;
  std::cout << "zm_eamxx_bridge_run - nlev: "<< nlev << std::endl;

  auto sz2 = T_mid.size() * Spack::n;

  std::cout << "zm_eamxx_bridge_run - sz1: "<< T_mid.size() << std::endl;
  std::cout << "zm_eamxx_bridge_run - sz2: "<< sz2 << std::endl;

  for (int i = 0; i<ncol; ++i) {
    // std::cout << "zm_eamxx_bridge_run - phis("<<i<<")     : " << phis_s(i)      << std::endl;
    std::cout << "zm_eamxx_bridge_run - T("<<i<<",     0): " << T_s(i,      0) << std::endl;
    std::cout << "zm_eamxx_bridge_run - q("<<i<<",     0): " << q_s(i,      0) << std::endl;
    std::cout << "zm_eamxx_bridge_run - T("<<i<<",  72-1): " << T_s(i,   72-1) << std::endl;
    std::cout << "zm_eamxx_bridge_run - q("<<i<<",  72-1): " << q_s(i,   72-1) << std::endl;
    std::cout << "zm_eamxx_bridge_run - T("<<i<<",nlev-1): " << T_s(i, nlev-1) << std::endl;
    std::cout << "zm_eamxx_bridge_run - q("<<i<<",nlev-1): " << q_s(i, nlev-1) << std::endl;
  }

  zm_input.transpose<ekat::TransposeDirection::c2f>();

  zm_eamxx_bridge_run_c( zm_input.ncol,
                         zm_input.is_first_step,
                         zm_input.phis.data(),
                         ekat::scalarize( zm_input.T_mid ).data(),
                         ekat::scalarize( zm_input.qv    ).data() );

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
