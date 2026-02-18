#ifndef ZM_INPUT_STATE_IMPL_HPP
#define ZM_INPUT_STATE_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

// -------------------------------------------------------------------------
// transpose method for fortran bridging
template<typename S, typename D>
template <ekat::TransposeDirection::Enum DirT>
void Functions<S,D>::ZmInputState::transpose(int ncol, int nlev_mid)
{
  auto nlev_int = nlev_mid+1;

  // ***********************************************************************
  // TEMPORARY
  // ***********************************************************************
  auto nlev_mid_packs = ekat::npack<Spack>(nlev_mid);
  auto nlev_int_packs = ekat::npack<Spack>(nlev_int);
  if (DirT == ekat::TransposeDirection::c2f) {
    // create temporaries to avoid "Implicit capture" warning
    const auto loc_f_z_mid   = f_z_mid;
    const auto loc_f_p_mid   = f_p_mid;
    const auto loc_f_p_del   = f_p_del;
    const auto loc_f_T_mid   = f_T_mid;
    const auto loc_f_qv      = f_qv;
    const auto loc_f_uwind   = f_uwind;
    const auto loc_f_vwind   = f_vwind;
    const auto loc_f_omega   = f_omega;
    const auto loc_f_cldfrac = f_cldfrac;
    const auto loc_f_z_int   = f_z_int;
    const auto loc_f_p_int   = f_p_int;

    //----------------------------------------------------------------------
    // mid-point level variables
    Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int i) {
      const int icol = i/nlev_mid_packs;
      const int klev = i%nlev_mid_packs;
      loc_f_z_mid   (icol,klev) = z_mid   (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_p_mid   (icol,klev) = p_mid   (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_p_del   (icol,klev) = p_del   (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_T_mid   (icol,klev) = T_mid   (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_qv      (icol,klev) = qv      (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_uwind   (icol,klev) = uwind   (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_vwind   (icol,klev) = vwind   (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_omega   (icol,klev) = omega   (icol,klev/Spack::n)[klev%Spack::n];
      loc_f_cldfrac (icol,klev) = cldfrac (icol,klev/Spack::n)[klev%Spack::n];
    });

    // interface level variables
    Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_int_packs), KOKKOS_LAMBDA (const int i) {
      const int icol = i/nlev_int_packs;
      const int klev = i%nlev_int_packs;
      f_z_int   (icol,klev) = z_int   (icol,klev/Spack::n)[klev%Spack::n];
      f_p_int   (icol,klev) = p_int   (icol,klev/Spack::n)[klev%Spack::n];
    });

    //----------------------------------------------------------------------
    // copy to host mirrors
    Kokkos::deep_copy(h_phis,     phis );
    Kokkos::deep_copy(h_pblh,     pblh );
    Kokkos::deep_copy(h_tpert,    tpert );
    Kokkos::deep_copy(h_landfrac, landfrac );
    Kokkos::deep_copy(h_z_mid,    f_z_mid );
    Kokkos::deep_copy(h_p_mid,    f_p_mid );
    Kokkos::deep_copy(h_p_del,    f_p_del );
    Kokkos::deep_copy(h_T_mid,    f_T_mid );
    Kokkos::deep_copy(h_qv,       f_qv );
    Kokkos::deep_copy(h_uwind,    f_uwind );
    Kokkos::deep_copy(h_vwind,    f_vwind );
    Kokkos::deep_copy(h_omega,    f_omega );
    Kokkos::deep_copy(h_cldfrac,  f_cldfrac );
    Kokkos::deep_copy(h_z_int,    f_z_int );
    Kokkos::deep_copy(h_p_int,    f_p_int );
  }

  // ***********************************************************************
  // TEMPORARY
  // ***********************************************************************

  // if (DirT == ekat::TransposeDirection::c2f) {
  //   ekat::device_to_host({h_phis.data()},     ncol,           std::vector< view_1d<const Scalar>>{phis});
  //   ekat::device_to_host({h_pblh.data()},     ncol,           std::vector< view_1d<const Scalar>>{pblh});
  //   ekat::device_to_host({h_tpert.data()},    ncol,           std::vector<uview_1d<      Scalar>>{tpert});
  //   ekat::device_to_host({h_landfrac.data()}, ncol,           std::vector< view_1d<const Scalar>>{landfrac});
  //   ekat::device_to_host({h_z_mid.data()},    ncol, nlev_mid, std::vector<uview_2d<const Spack >>{z_mid});
  //   ekat::device_to_host({h_p_mid.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{p_mid});
  //   ekat::device_to_host({h_p_del.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{p_del});
  //   ekat::device_to_host({h_T_mid.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{T_mid});
  //   ekat::device_to_host({h_qv.data()},       ncol, nlev_mid, std::vector< view_2d<      Spack >>{qv});
  //   ekat::device_to_host({h_uwind.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{uwind});
  //   ekat::device_to_host({h_vwind.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{vwind});
  //   ekat::device_to_host({h_omega.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{omega});
  //   ekat::device_to_host({h_cldfrac.data()},  ncol, nlev_mid, std::vector< view_2d<const Spack >>{cldfrac});
  //   ekat::device_to_host({h_z_int.data()},    ncol, nlev_int, std::vector<uview_2d<      Spack >>{z_int});
  //   ekat::device_to_host({h_p_int.data()},    ncol, nlev_int, std::vector< view_2d<const Spack >>{p_int});
  // }
}

template<typename S, typename D>
void Functions<S,D>::ZmInputState::calculate_tpert(int ncol, int nlev, bool is_first_step)
{
  const Real cpair  = PC::Cpair.value;
  const Real latvap = PC::LatVap.value;

  // create temporaries to avoid "Implicit capture" warning
  auto loc_tpert    = tpert;
  auto loc_pblh     = pblh;
  auto loc_z_int    = z_int;
  auto loc_p_mid    = p_mid;
  auto loc_qc       = qc;
  auto loc_thl_sec  = thl_sec;

  Kokkos::parallel_for("zm_calculate_tpert",ncol, KOKKOS_LAMBDA (const int i) {
    if (is_first_step) {
      loc_tpert(i) = 0.0;
    } else {
      // identify interface index for top of PBL
      int pblh_k_ind = -1;
      for (int k=0; k<nlev; ++k) {
        auto z_int_tmp_k   = loc_z_int(i,k/Spack::n)[k%Spack::n];
        auto z_int_tmp_kp1 = loc_z_int(i,k/Spack::n)[k%Spack::n];
        if ( z_int_tmp_k>loc_pblh(i) && z_int_tmp_kp1<=loc_pblh(i) ) {
          pblh_k_ind = k;
        }
      }
      if (pblh_k_ind==-1) {
        // PBL top index not found, so just set the perturbation to zero
        loc_tpert(i) = 0.0;
      } else {
        // calculate tpert as std deviation of temperature from SHOC's theta-l variance
        auto exner_pbl    = PF::exner_function( loc_p_mid(i,pblh_k_ind/Spack::n)[pblh_k_ind%Spack::n] );
        auto qc_pbl       = loc_qc(i,pblh_k_ind/Spack::n)[pblh_k_ind%Spack::n];
        auto thl_sec_pbl  = loc_thl_sec(i,pblh_k_ind/Spack::n)[pblh_k_ind%Spack::n];
        auto thl_std_pbl  = sqrt( thl_sec_pbl ); // std deviation of thetal;
        loc_tpert(i) = ( thl_std_pbl + (latvap/cpair)*qc_pbl ) / exner_pbl;
        loc_tpert(i) = ekat::impl::min(2.0,loc_tpert(i)); // apply limiter
      }
    }
  });
}

} // namespace zm
} // namespace scream

#endif
