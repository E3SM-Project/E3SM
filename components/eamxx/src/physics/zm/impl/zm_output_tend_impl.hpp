#ifndef ZM_OUTPUT_TEND_IMPL_HPP
#define ZM_OUTPUT_TEND_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

// -------------------------------------------------------------------------
// transpose method for fortran bridging
template<typename S, typename D>
template <ekat::TransposeDirection::Enum DirT>
void Functions<S,D>::ZmOutputTend::transpose(int ncol, int nlev_mid)
{
  auto nlev_int = nlev_mid+1;

  // ***********************************************************************
  // TEMPORARY
  // ***********************************************************************
  auto nlev_mid_packs = ekat::npack<Spack>(nlev_mid);
  auto nlev_int_packs = ekat::npack<Spack>(nlev_int);
  if (DirT == ekat::TransposeDirection::f2c) {
    // copy back to device
    Kokkos::deep_copy(f_tend_t,   h_tend_t);
    Kokkos::deep_copy(f_tend_qv,  h_tend_qv);
    Kokkos::deep_copy(f_tend_u,   h_tend_u);
    Kokkos::deep_copy(f_tend_v,   h_tend_v);
    Kokkos::deep_copy(f_rain_prod,h_rain_prod);
    Kokkos::deep_copy(f_snow_prod,h_snow_prod);
    Kokkos::deep_copy(f_prec_flux,h_prec_flux);
    Kokkos::deep_copy(f_snow_flux,h_snow_flux);
    Kokkos::deep_copy(f_mass_flux,h_mass_flux);
    Kokkos::deep_copy(prec,       h_prec);
    Kokkos::deep_copy(snow,       h_snow);
    Kokkos::deep_copy(cape,       h_cape);
    Kokkos::deep_copy(activity,   h_activity);

    //----------------------------------------------------------------------
    // create temporaries to avoid "Implicit capture" warning
    const auto loc_tend_t    = tend_t;
    const auto loc_tend_qv   = tend_qv;
    const auto loc_tend_u    = tend_u;
    const auto loc_tend_v    = tend_v;
    const auto loc_rain_prod = rain_prod;
    const auto loc_snow_prod = snow_prod;
    const auto loc_prec_flux = prec_flux;
    const auto loc_snow_flux = snow_flux;
    const auto loc_mass_flux = mass_flux;

    //----------------------------------------------------------------------
    // mid-point level variables
    Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int i) {
      const int icol = i/nlev_mid_packs;
      const int klev = i%nlev_mid_packs;
      loc_tend_t   (icol,klev/Spack::n)[klev%Spack::n] = f_tend_t   (icol,klev);
      loc_tend_qv  (icol,klev/Spack::n)[klev%Spack::n] = f_tend_qv  (icol,klev);
      loc_tend_u   (icol,klev/Spack::n)[klev%Spack::n] = f_tend_u   (icol,klev);
      loc_tend_v   (icol,klev/Spack::n)[klev%Spack::n] = f_tend_v   (icol,klev);
      loc_rain_prod(icol,klev/Spack::n)[klev%Spack::n] = f_rain_prod(icol,klev);
      loc_snow_prod(icol,klev/Spack::n)[klev%Spack::n] = f_snow_prod(icol,klev);
    });

    // interface level variables
    Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_int_packs), KOKKOS_LAMBDA (const int i) {
      const int icol = i/nlev_int_packs;
      const int klev = i%nlev_int_packs;
      loc_prec_flux(icol,klev/Spack::n)[klev%Spack::n] = f_prec_flux(icol,klev);
      loc_snow_flux(icol,klev/Spack::n)[klev%Spack::n] = f_snow_flux(icol,klev);
      loc_mass_flux(icol,klev/Spack::n)[klev%Spack::n] = f_mass_flux(icol,klev);
    });
  }
  // ***********************************************************************
  // TEMPORARY
  // ***********************************************************************

  // if (DirT == ekat::TransposeDirection::f2c) {
  //   ekat::host_to_device({h_prec.data()},     ncol,           std::vector<uview_1d<Scalar>>{prec});
  //   ekat::host_to_device({h_activity.data()}, ncol,           std::vector<uview_1d<Int>>   {activity});
  //   ekat::host_to_device({h_prec.data()},     ncol,           std::vector<uview_1d<Scalar>>{prec});
  //   ekat::host_to_device({h_snow.data()},     ncol,           std::vector<uview_1d<Scalar>>{snow});
  //   ekat::host_to_device({h_cape.data()},     ncol,           std::vector<uview_1d<Scalar>>{cape});
  //   ekat::host_to_device({h_tend_t.data()},   ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_t},    true);
  //   ekat::host_to_device({h_tend_qv.data()},  ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_qv},   true);
  //   ekat::host_to_device({h_tend_u.data()},   ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_u},    true);
  //   ekat::host_to_device({h_tend_v.data()},   ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_v},    true);
  //   ekat::host_to_device({h_rain_prod.data()},ncol, nlev_mid, std::vector<uview_2d<Spack>> {rain_prod}, true);
  //   ekat::host_to_device({h_snow_prod.data()},ncol, nlev_mid, std::vector<uview_2d<Spack>> {snow_prod}, true);
  //   ekat::host_to_device({h_prec_flux.data()},ncol, nlev_int, std::vector<uview_2d<Spack>> {prec_flux}, true);
  //   ekat::host_to_device({h_snow_flux.data()},ncol, nlev_int, std::vector<uview_2d<Spack>> {snow_flux}, true);
  //   ekat::host_to_device({h_mass_flux.data()},ncol, nlev_int, std::vector<uview_2d<Spack>> {mass_flux}, true);
  // }
}

template<typename S, typename D>
void Functions<S,D>::ZmOutputTend::init(int ncol, int nlev_mid)
{
  auto nlev_int = nlev_mid+1;
  auto nlev_mid_packs = ekat::npack<Spack>(nlev_mid);
  auto nlev_int_packs = ekat::npack<Spack>(nlev_int);
  Real init_fill_value = 0;
  // create temporaries to avoid "Implicit capture" warning
  auto loc_prec       = prec;
  auto loc_snow       = snow;
  auto loc_cape       = cape;
  auto loc_activity   = activity;
  auto loc_tend_t     = tend_t;
  auto loc_tend_qv    = tend_qv;
  auto loc_tend_u     = tend_u;
  auto loc_tend_v     = tend_v;
  auto loc_rain_prod  = rain_prod;
  auto loc_snow_prod  = snow_prod;
  auto loc_prec_flux  = prec_flux;
  auto loc_snow_flux  = snow_flux;
  auto loc_mass_flux  = mass_flux;

  // 1D scalar variables
  Kokkos::parallel_for("zm_output_init_s", KT::RangePolicy(0, ncol), KOKKOS_LAMBDA (const int i) {
    loc_prec(i)     = init_fill_value;
    loc_snow(i)     = init_fill_value;
    loc_cape(i)     = init_fill_value;
    loc_activity(i) = -1;
  });

  // mid-point level variables
  Kokkos::parallel_for("zm_output_init_m",KT::RangePolicy(0, ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int i) {
    const int icol = i/nlev_mid_packs;
    const int klev = i%nlev_mid_packs;
    loc_tend_t   (icol,klev) = init_fill_value;
    loc_tend_qv  (icol,klev) = init_fill_value;
    loc_tend_u   (icol,klev) = init_fill_value;
    loc_tend_v   (icol,klev) = init_fill_value;
    loc_rain_prod(icol,klev) = init_fill_value;
    loc_snow_prod(icol,klev) = init_fill_value;
  });

  // interface level variables
  Kokkos::parallel_for("zm_output_init_i",KT::RangePolicy(0, ncol*nlev_int_packs), KOKKOS_LAMBDA (const int i) {
    const int icol = i/nlev_int_packs;
    const int klev = i%nlev_int_packs;
    loc_prec_flux(icol,klev) = init_fill_value;
    loc_snow_flux(icol,klev) = init_fill_value;
    loc_mass_flux(icol,klev) = init_fill_value;
  });
}

} // namespace zm
} // namespace scream

#endif
