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
                        const Real *state_p_mid,
                        const Real *state_p_int,
                        const Real *state_p_del,
                              Real *state_t,
                              Real *state_qv,
                              Real *state_qc,
                              Real *state_omega,
                        const Real *state_pblh,
                              Real *output_prec,
                              Real *output_cape,
                              Real *output_tend_s,
                              Real *output_tend_q,
                              Real *output_prec_flux,
                              Real *output_mass_flux
                            );
} // extern "C" : end _c decls

namespace scream {
namespace zm {

void zm_eamxx_bridge_init( Int pcol, Int pver ){
  zm_eamxx_bridge_init_c( pcol, pver );
}

void zm_eamxx_bridge_run( Int pver,
                          ZMF::zm_input_state& zm_input,
                          ZMF::zm_output_tend& zm_output
                        ){
  //----------------------------------------------------------------------------
  // auto s_p = ekat::scalarize( zm_input.p_mid );
  // auto s_T = ekat::scalarize( zm_input.T_mid );
  // auto s_q = ekat::scalarize( zm_input.qv    );
  // for (int i=0; i<zm_input.ncol; ++i) {
  //   for (int k=0; k<12; ++k) {
  //     std::cout << "zm_eamxx_bridge_run - p("<<i<<","<<k<<"): " << s_p(i,k) << std::endl;
  //   }
  // }
  //----------------------------------------------------------------------------
  zm_input.transpose<ekat::TransposeDirection::c2f>(pver);
  zm_output.transpose<ekat::TransposeDirection::c2f>(pver);

  for (int i=0; i<zm_output.ncol; ++i) {
    for (int k=0; k<pver; ++k) {
      std::cout << "zm_eamxx_bridge_run (pre fortran) - ("<<i<<","<<k<<") - "
      << "precip "    << " / "
      <<"cape "       << " / "
      <<"f_tend_s "   << " / "
      <<"f_tend_q "
      <<" : "
      << zm_output.precip(i)      << " / "
      << zm_output.cape(i)        << " / "
      << zm_output.f_tend_s(i,k)  << " / "
      << zm_output.f_tend_q(i,k)
      << std::endl;
    }
  }

  zm_eamxx_bridge_run_c( zm_input.ncol,
                         zm_input.dtime,
                         zm_input.is_first_step,
                         zm_input.phis          .data(),
                         zm_input.f_p_mid       .data(),
                         zm_input.f_p_int       .data(),
                         zm_input.f_p_del       .data(),
                         zm_input.f_T_mid       .data(),
                         zm_input.f_qv          .data(),
                         zm_input.f_qc          .data(),
                         zm_input.f_omega       .data(),
                         zm_input.pblh          .data(),
                         zm_output.precip       .data(),
                         zm_output.cape         .data(),
                         zm_output.f_tend_s     .data(),
                         zm_output.f_tend_q     .data(),
                         zm_output.f_prec_flux  .data(),
                         zm_output.f_mass_flux  .data()
                        );

  for (int i=0; i<zm_output.ncol; ++i) {
    // std::cout << "zm_eamxx_bridge_run (post fortran) - ("<<i<<") - "<<"precip   : "<< zm_output.precip(i)     << std::endl;
    // std::cout << "zm_eamxx_bridge_run (post fortran) - ("<<i<<") - "<<"cape     : "<< zm_output.cape(i)       << std::endl;
    for (int k=0; k<pver; ++k) {
      // std::cout << "zm_eamxx_bridge_run (post fortran) - ("<<i<<","<<k<<") - "<< std::endl;

      // std::cout << "zm_eamxx_bridge_run (post fortran) - ("<<i<<","<<k<<") - "
      // <<"precip "     << " / "
      // <<"cape "       << " / "
      // <<"f_tend_s "   << " / "
      // <<"f_tend_q "
      // <<" : "
      // << zm_output.precip(i)      << " / "
      // << zm_output.cape(i)        << " / "
      // << zm_output.f_tend_s(i,k)  << " / "
      // << zm_output.f_tend_q(i,k)
      // << std::endl;

      std::cout << "zm_eamxx_bridge_run (post fortran) - ("<<i<<","<<k<<") - "<<"f_tend_s : "<< zm_output.f_tend_s(i,k) << std::endl;
      std::cout << "zm_eamxx_bridge_run (post fortran) - ("<<i<<","<<k<<") - "<<"f_tend_q : "<< zm_output.f_tend_q(i,k) << std::endl;
    }
  }

  zm_input.transpose<ekat::TransposeDirection::f2c>(pver);
  zm_output.transpose<ekat::TransposeDirection::f2c>(pver);

  //----------------------------------------------------------------------------
  // auto s_p_mid  = ekat::scalarize( zm_input.p_mid );
  // auto s_tend_s = ekat::scalarize( zm_output.tend_s );
  // auto s_tend_q = ekat::scalarize( zm_output.tend_q );

  // auto pcol = zm_output.pcol;

  // auto tend_s_tmp = view_2d<Real>("tend_s_tmp", pcol, pver);
  // auto tend_q_tmp = view_2d<Real>("tend_q_tmp", pcol, pver);
  // for (int i=0; i<pcol; ++i) {
  //   for (int j=0; j<pver; ++j) {
  //     tend_s_tmp(i,j) = zm_outputtend_s(i, j / Spack::n)[j % Spack::n];
  //     tend_q_tmp(i,j) = zm_outputtend_q(i, j / Spack::n)[j % Spack::n];
  //   }
  // }

  // for (int i=0; i<zm_input.ncol; ++i) {
  //   for (int k=0; k<pver; ++k) {
  //     std::cout << "zm_eamxx_bridge_run - ("<<i<<","<<k<<") f_tend_s / tend_s / f_tend_q / tend_q: " 
  //     << zm_output.f_tend_s(i,k) << " / " << zm_output.  tend_s(i,k/16)[k%16] << " / " 
  //     << zm_output.f_tend_q(i,k) << " / " << zm_output.  tend_q(i,k/16)[k%16] << " / " 
  //     << std::endl;
  //   }
  // }

  // for (int i=0; i<zm_input.ncol; ++i) {
    
    // for (int k=0; k<6; ++k) {
    //   // std::cout << "zm_eamxx_bridge_run - ("<<i<<","<<k<<") tend_s / tend_q : " << s_tend_s(i,k) << " / " << s_tend_q(i,k) << std::endl;
    //   std::cout << "zm_eamxx_bridge_run - ("<<i<<","<<k<<") tend_s / tend_q : " << zm_output.tend_s(i,k) << " / " << zm_output.tend_q(i,k) << std::endl;
    // }
  // }

  // std::cout << " s_tend_s.extent(0): " << s_tend_s.extent(0) <<std::endl;
  // std::cout << " s_tend_q.extent(1): " << s_tend_q.extent(1) <<std::endl;

  // for (int i=0; i<s_tend_s.extent(0); ++i) {
  //   // std::cout << "zm_eamxx_bridge_run - prec("<<i<<"): " << zm_output.precip(i) << std::endl;
  //   for (int k=0; k<6; ++k) {
  //     std::cout << "zm_eamxx_bridge_run - ("<<i<<","<<k<<")  tend_s / tend_q : " << tend_s(i,k) << " / " << tend_q(i,k) << std::endl;
  //   }
  // }
  //----------------------------------------------------------------------------
}

// end _c impls

} // namespace zm
} // namespace scream
