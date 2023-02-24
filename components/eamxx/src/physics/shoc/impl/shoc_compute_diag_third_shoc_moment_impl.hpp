#ifndef SHOC_COMPUTE_DIAG_THIRD_SHOC_MOMENT_IMPL_HPP
#define SHOC_COMPUTE_DIAG_THIRD_SHOC_MOMENT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::compute_diag_third_shoc_moment(
  const MemberType& team,
  const Int& nlev,
  const Int& nlevi,
  const uview_1d<const Spack>& w_sec,
  const uview_1d<const Spack>& thl_sec,
  const uview_1d<const Spack>& wthl_sec,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& isotropy_zi,
  const uview_1d<const Spack>& brunt_zi,
  const uview_1d<const Spack>& w_sec_zi,
  const uview_1d<const Spack>& thetal_zi,
  const uview_1d<Spack>& w3)
{
  const auto ggr = C::gravit;

  const Int nlev_pack = ekat::npack<Spack>(nlev);

  // Scalarize views for shifts and single entry access
  const auto s_dz_zt = scalarize(dz_zt);
  const auto s_wthl_sec = scalarize(wthl_sec);
  const auto s_thl_sec = scalarize(thl_sec);
  const auto s_w_sec = scalarize(w_sec);
  const auto s_tke = scalarize(tke);
  const auto s_w3 = scalarize(w3);

  const Int max_safe_shift1_idx = nlevi - 2;

  // Set lower condition: w3(i,nlevi) = 0
  s_w3(nlevi-1) = 0;

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    // Constants
    const auto c_diag_3rd_mom = scream::shoc::Constants<Scalar>::c_diag_3rd_mom;
    const Scalar a0 = (sp(0.52)*(1/(c_diag_3rd_mom*c_diag_3rd_mom)))/(c_diag_3rd_mom-2);
    const Scalar a1 = sp(0.87)/(c_diag_3rd_mom*c_diag_3rd_mom);
    const Scalar a2 = sp(0.5)/c_diag_3rd_mom;
    const Scalar a3 = sp(0.6)/(c_diag_3rd_mom*(c_diag_3rd_mom-2));
    const Scalar a4 = sp(2.4)/(3*c_diag_3rd_mom+5);
    const Scalar a5 = sp(0.6)/(c_diag_3rd_mom*(3+5*c_diag_3rd_mom));

    // Calculate shifts
    Spack
      dz_zt_k,    dz_zt_km1,
      wthl_sec_k, wthl_sec_km1, wthl_sec_kp1,
      thl_sec_k,  thl_sec_km1,  thl_sec_kp1,
      w_sec_k,    w_sec_km1,
      tke_k,      tke_km1;

    auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    auto range_pack_m1 = range_pack;
    auto range_pack_m2 = range_pack;
    // index for _km1 should never go below 0
    range_pack_m1.set(range_pack < 1, 1);
    // index for _km1 should never go above scalar view bounds
    range_pack_m2.set(range_pack > max_safe_shift1_idx, max_safe_shift1_idx);

    ekat::index_and_shift<-1>(s_dz_zt, range_pack_m1, dz_zt_k, dz_zt_km1);
    ekat::index_and_shift<-1>(s_wthl_sec, range_pack_m1, wthl_sec_k, wthl_sec_km1);
    ekat::index_and_shift<1> (s_wthl_sec, range_pack_m2, wthl_sec_k, wthl_sec_kp1);
    ekat::index_and_shift<-1>(s_thl_sec, range_pack_m1, thl_sec_k, thl_sec_km1);
    ekat::index_and_shift<1> (s_thl_sec, range_pack_m2, thl_sec_k, thl_sec_kp1);
    ekat::index_and_shift<-1>(s_w_sec, range_pack_m1, w_sec_k, w_sec_km1);
    ekat::index_and_shift<-1>(s_tke, range_pack_m1, tke_k, tke_km1);

    const auto active_range = range_pack > 0 && range_pack < nlev;
    if (active_range.any()) {
      // Compute inputs for computing f0 to f5 terms
      const auto thedz  = 1/dz_zi(k);
      const auto thedz2 = 1/(dz_zt_k+dz_zt_km1);

      const auto iso       = isotropy_zi(k);
      const auto isosqrd   = ekat::square(iso);
      const auto buoy_sgs2 = isosqrd*brunt_zi(k);
      const auto bet2      = ggr/thetal_zi(k);

      // Compute f0 to f5 terms
      const Spack thl_sec_diff = thl_sec_km1 - thl_sec_kp1;
      const Spack wthl_sec_diff = wthl_sec_km1 - wthl_sec_kp1;
      const Spack wsec_diff = w_sec_km1 - w_sec(k);
      const Spack tke_diff = tke_km1 - tke(k);

      const auto f0 = thedz2*ekat::cube(bet2)*((iso*iso)*(iso*iso))*wthl_sec_k*thl_sec_diff;
      const auto f1 = thedz2*ekat::square(bet2)*ekat::cube(iso)*(wthl_sec_k*wthl_sec_diff+sp(0.5)*w_sec_zi(k)*thl_sec_diff);
      const auto f2 = thedz*bet2*isosqrd*wthl_sec_k*wsec_diff+2*thedz2*bet2*isosqrd*w_sec_zi(k)*wthl_sec_diff;
      const auto f3 = thedz2*bet2*isosqrd*w_sec_zi(k)*wthl_sec_diff+thedz*bet2*isosqrd*(wthl_sec_k*tke_diff);
      const auto f4 = thedz*iso*w_sec_zi(k)*(wsec_diff+tke_diff);
      const auto f5 = thedz*iso*w_sec_zi(k)*wsec_diff;

      // Compute omega terms
      const auto omega0 = a4/Spack(1-a5*buoy_sgs2);
      const auto omega1 = omega0/(2*c_diag_3rd_mom);
      const auto omega2 = omega1*f3+sp(5.0/4.0)*omega0*f4;

      // Compute the x0, y0, x1, y1 terms
      const auto x0 = (a2*buoy_sgs2*(Spack(1)-a3*buoy_sgs2))/(Spack(1)-(a1+a3)*buoy_sgs2);
      const auto y0 = (2*a2*buoy_sgs2*x0)/(Spack(1)-a3*buoy_sgs2);
      const auto x1 = (a0*f0+a1*f1+a2*(Spack(1)-a3*buoy_sgs2)*f2)/(Spack(1)-(a1+a3)*buoy_sgs2);
      const auto y1 = (2*a2*(buoy_sgs2*x1+(a0/a1)*f0+f1))/(Spack(1)-a3*buoy_sgs2);

      // Compute the aa0, aa1 terms
      const auto aa0 = omega0*x0+omega1*y0;
      const auto aa1 = omega0*x1+omega1*y1+omega2;

      // Finally, compute the third moment of w
      w3(k).set(active_range,
                (aa1-sp(1.2)*x1-sp(1.5)*f5)/(Spack(c_diag_3rd_mom)-sp(1.2)*x0+aa0));
    }
  });

  // Set upper condition: w3(i,0) = 0
  s_w3(0) = 0;
}

} // namespace shoc
} // namespace scream

#endif
