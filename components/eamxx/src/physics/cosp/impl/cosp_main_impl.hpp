#ifndef SCREAM_COSP_MAIN_IMPL_HPP
#define SCREAM_COSP_MAIN_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

#include <ekat_assert.hpp>

namespace scream {
namespace cosp {

/*
 * The C++ counterpart of cosp_c2f_run (cosp_c2f.F90): the optical inputs of
 * subsample_and_optics, the input checks of cosp_errorCheck that matter for the
 * ISCCP/MISR/MODIS outputs, the subcolumn and column simulators of cosp_simulator, and
 * the reordering of the outputs. Everything stays on device.
 */

template <typename S, typename D>
Functions<S,D>::Workspace::Workspace(const int ncol, const int nsubcol, const int nlev)
 : tau        ("cosp_ws_tau",         ncol, nsubcol, nlev)
 , dem        ("cosp_ws_dem",         ncol, nsubcol, nlev)
 , liqfrac    ("cosp_ws_liqfrac",     ncol, nsubcol, nlev)
 , g          ("cosp_ws_g",           ncol, nsubcol, nlev)
 , w0         ("cosp_ws_w0",          ncol, nsubcol, nlev)
 , boxtau     ("cosp_ws_boxtau",      ncol, nsubcol)
 , boxptop    ("cosp_ws_boxptop",     ncol, nsubcol)
 , boxttop    ("cosp_ws_boxttop",     ncol, nsubcol)
 , misr_tau   ("cosp_ws_misr_tau",    ncol, nsubcol)
 , misr_ztop  ("cosp_ws_misr_ztop",   ncol, nsubcol)
 , modis_ctp  ("cosp_ws_modis_ctp",   ncol, nsubcol)
 , modis_tau  ("cosp_ws_modis_tau",   ncol, nsubcol)
 , modis_size ("cosp_ws_modis_size",  ncol, nsubcol)
 , modis_phase("cosp_ws_modis_phase", ncol, nsubcol)
{}

template <typename S, typename D>
void Functions<S,D>::cosp_main(
    const int ncol, const int nsubcol, const int nlev, const S emsfc_lw,
    const view_1d<const int>& sunlit, const view_1d<const S>& skt,
    const view_2d<const S>& T_mid, const view_2d<const S>& p_mid,
    const view_2d<const S>& p_int, const view_2d<const S>& z_mid,
    const view_2d<const S>& qv, const view_2d<const S>& qc, const view_2d<const S>& qi,
    const view_2d<const S>& cldfrac, const view_2d<const S>& reff_qc,
    const view_2d<const S>& reff_qi, const view_2d<const S>& dtau067,
    const view_2d<const S>& dtau105, const view_1d<S>& isccp_cldtot,
    const view_3d<S>& isccp_ctptau, const view_3d<S>& modis_ctptau,
    const view_3d<S>& misr_cthtau, const Workspace& ws)
{
  EKAT_REQUIRE_MSG(nsubcol == 1,
      "Error! The C++ COSP only supports cosp_subcolumns=1 so far (got " +
      std::to_string(nsubcol) + ").\n");

  using ExeSpace = typename KT::ExeSpace;
  using Policy2  = Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<2>>;
  using Policy3  = Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>>;

  constexpr int ntau   = C::ntau;
  constexpr int npres  = C::npres;
  constexpr int nhgt   = C::nhgt;
  constexpr S   undef  = C::R_UNDEF;

  // Bits of the simulators that cosp_errorCheck switches off
  constexpr int isccp_off = 1, misr_off = 2, modis_off = 4;

  const auto tau = ws.tau, dem = ws.dem, liqfrac = ws.liqfrac, g = ws.g, w0 = ws.w0;
  const auto boxtau = ws.boxtau, boxptop = ws.boxptop, boxttop = ws.boxttop;
  const auto misr_tau = ws.misr_tau, misr_ztop = ws.misr_ztop;
  const auto modis_ctp = ws.modis_ctp, modis_tau = ws.modis_tau, modis_size = ws.modis_size;
  const auto modis_phase = ws.modis_phase;

  // -----------------------------------------------------------------------------------
  // Optical inputs (subsample_and_optics). With one subcolumn the inputs are converted
  // from in-cloud to grid-box values, and frac_out is 1 (stratiform) everywhere.
  Kokkos::parallel_for("cosp_optics", Policy3({0,0,0},{ncol,nsubcol,nlev}),
                       KOKKOS_LAMBDA (const int j, const int i, const int k) {
    const S frac_out = C::I_LSC;
    const S dtau_s = cldfrac(j,k)*dtau067(j,k);
    const S dem_s  = 1 - Kokkos::exp(-cldfrac(j,k)*dtau105(j,k));
    const S t = simulator_optics(frac_out, 0, dtau_s);
    tau(j,i,k) = t;
    dem(j,i,k) = simulator_optics(frac_out, 0, dem_s);

    // MODIS: there is no convective cloud, and the sizes go through m and back
    const S cloud_water = simulator_optics(frac_out, 0, qc(j,k));
    const S cloud_ice   = simulator_optics(frac_out, 0, qi(j,k));
    const S water_size  = simulator_optics(frac_out, 0, S(1e-6)*reff_qc(j,k));
    const S ice_size    = simulator_optics(frac_out, 0, S(1e-6)*reff_qi(j,k));
    S tau_liq, tau_ice;
    modis_optics_partition(cloud_water, cloud_ice, water_size, ice_size, t, tau_liq, tau_ice);
    modis_optics(tau_liq, water_size*S(1.0e6), tau_ice, ice_size*S(1.0e6),
                 liqfrac(j,i,k), g(j,i,k), w0(j,i,k));
  });

  // -----------------------------------------------------------------------------------
  // Input checks (cosp_errorCheck): an out-of-range value anywhere switches a simulator
  // off for all columns, and its outputs are set to R_UNDEF
  int off = 0;
  Kokkos::parallel_reduce("cosp_check_state", Policy2({0,0},{ncol,nlev+1}),
                          KOKKOS_LAMBDA (const int j, const int k, int& f) {
    if (k == 0) {
      if (sunlit(j) < 0) f |= isccp_off | misr_off | modis_off;
      if (skt(j) < 0)    f |= isccp_off;
    }
    if (p_int(j,k) < 0) f |= isccp_off | modis_off;
    if (k == nlev) return;
    if (T_mid(j,k) < 0)    f |= isccp_off | misr_off;
    if (p_mid(j,k) < 0)    f |= isccp_off;
    if (qv(j,k) < 0)       f |= isccp_off;
    if (z_mid(j,k) < -300) f |= misr_off;
  }, Kokkos::BOr<int>(off));
  int off_optics = 0;
  Kokkos::parallel_reduce("cosp_check_optics", Policy3({0,0,0},{ncol,nsubcol,nlev}),
                          KOKKOS_LAMBDA (const int j, const int i, const int k, int& f) {
    if (tau(j,i,k) < 0)                     f |= isccp_off | misr_off | modis_off;
    if (dem(j,i,k) < 0 || dem(j,i,k) > 1)   f |= isccp_off;
    if (g(j,i,k) < -1 || g(j,i,k) > 1)      f |= modis_off;
    if (w0(j,i,k) < 0 || w0(j,i,k) > 1)     f |= modis_off;
  }, Kokkos::BOr<int>(off_optics));
  off |= off_optics;
  if (emsfc_lw < 0 || emsfc_lw > 1) off |= isccp_off;

  const bool do_isccp = !(off & isccp_off);
  const bool do_misr  = !(off & misr_off);
  const bool do_modis = !(off & modis_off);

  // -----------------------------------------------------------------------------------
  // Simulators, one column per thread
  Kokkos::parallel_for("cosp_simulators", Kokkos::RangePolicy<ExeSpace>(0,ncol),
                       KOKKOS_LAMBDA (const int j) {
    using Kokkos::ALL;
    using Kokkos::subview;
    const int sl = sunlit(j);
    const auto qv_j = subview(qv,    j, ALL);
    const auto at_j = subview(T_mid, j, ALL);
    const auto pf_j = subview(p_mid, j, ALL);
    const auto ph_j = subview(p_int, j, ALL);
    const auto zf_j = subview(z_mid, j, ALL);

    // ISCCP subcolumns, which MODIS also needs
    if (do_isccp || do_modis) {
      const auto pt = icarus_point(nlev, skt(j), emsfc_lw, qv_j, at_j, pf_j, ph_j);
      for (int i = 0; i < nsubcol; ++i) {
        icarus_subcolumn(nlev, sl, pt, subview(tau, j, i, ALL), subview(dem, j, i, ALL),
                         skt(j), emsfc_lw, qv_j, at_j, pf_j, ph_j,
                         boxtau(j,i), boxptop(j,i), boxttop(j,i));
      }
    }

    // ISCCP column: the pressure bins are reordered to go from the surface up
    const auto isccp_fq = subview(isccp_ctptau, j, ALL, ALL);
    if (do_isccp) {
      S cldarea;
      icarus_column(nsubcol, sl, subview(boxtau, j, ALL), subview(boxptop, j, ALL),
                    isccp_fq, cldarea);
      for (int it = 0; it < ntau; ++it) {
        for (int ip = 0; ip < npres/2; ++ip) {
          const S tmp = isccp_fq(it,ip);
          isccp_fq(it,ip) = isccp_fq(it,npres-1-ip);
          isccp_fq(it,npres-1-ip) = tmp;
        }
      }
      // Values slightly off 1; the bounds are single precision in the Fortran
      if (cldarea > S(1.0f-1.e-5f) && cldarea < S(1.0f+1.e-5f)) cldarea = 1;
      isccp_cldtot(j) = cldarea;
    } else {
      for (int it = 0; it < ntau; ++it)
        for (int ip = 0; ip < npres; ++ip)
          isccp_fq(it,ip) = undef;
      isccp_cldtot(j) = undef;
    }

    // MISR
    const auto misr_fq = subview(misr_cthtau, j, ALL, ALL);
    if (do_misr) {
      for (int i = 0; i < nsubcol; ++i) {
        misr_subcolumn(nlev, subview(tau, j, i, ALL), zf_j, at_j, sl, misr_tau(j,i),
                       misr_ztop(j,i));
      }
      misr_column(nsubcol, subview(misr_ztop, j, ALL), sl, subview(misr_tau, j, ALL), misr_fq);
    } else {
      for (int it = 0; it < ntau; ++it)
        for (int ih = 0; ih < nhgt; ++ih)
          misr_fq(it,ih) = undef;
    }

    // MODIS, sunlit points only; the pressure bins are reordered as for ISCCP
    const auto modis_fq = subview(modis_ctptau, j, ALL, ALL);
    if (do_modis && sl > 0) {
      for (int i = 0; i < nsubcol; ++i) {
        modis_subcolumn(nlev, ph_j, subview(tau, j, i, ALL), subview(liqfrac, j, i, ALL),
                        subview(g, j, i, ALL), subview(w0, j, i, ALL), boxptop(j,i),
                        modis_phase(j,i), modis_ctp(j,i), modis_tau(j,i), modis_size(j,i));
      }
      modis_column(nsubcol, subview(modis_phase, j, ALL), subview(modis_ctp, j, ALL),
                   subview(modis_tau, j, ALL), subview(modis_size, j, ALL), modis_fq);
      for (int it = 0; it < ntau; ++it) {
        for (int ip = 0; ip < npres/2; ++ip) {
          const S tmp = modis_fq(it,ip);
          modis_fq(it,ip) = modis_fq(it,npres-1-ip);
          modis_fq(it,npres-1-ip) = tmp;
        }
      }
    } else {
      for (int it = 0; it < ntau; ++it)
        for (int ip = 0; ip < npres; ++ip)
          modis_fq(it,ip) = undef;
    }
  });
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_MAIN_IMPL_HPP
