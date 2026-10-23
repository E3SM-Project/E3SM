#include <catch2/catch.hpp>

#include "physics/cosp/tests/cosp_test_utils.hpp"

#include <cmath>

/*
 * BFB tests of the C++ MODIS simulator against the Fortran COSP library.
 */

using namespace scream;
using namespace scream::cosp_test;

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_two_stream_bfb", "[cosp]")
{
  // Random layers, plus single scattering albedos right at the (single precision)
  // conservative scattering limit, and very thick layers that hit the exp cutoff
  const int n = 2000;
  lview_1d tau("tau", n), g("g", n), w0("w0", n);
  fill_random(tau, 0, 50);
  fill_random(g, 0.6, 0.95);
  fill_random(w0, 0.9, 1.0);
  const double w0_cons = static_cast<double>(0.9999999f);
  const double edge_w0[] = {w0_cons, std::nextafter(w0_cons, 0.0), std::nextafter(w0_cons, 1.0),
                            1.0};
  for (int i = 0; i < 4; ++i) w0(i) = edge_w0[i];
  for (int i = 4; i < 20; ++i) tau(i) = uniform(1e3, 1e5);

  lview_1d ref_f("ref_f", n), tra_f("tra_f", n), ref2_f("ref2_f", n);
  cosp_two_stream_f(n, tau.data(), g.data(), w0.data(), ref_f.data(), tra_f.data(),
                    ref2_f.data());

  auto tau_d = to_device<dview_1d>(tau, "tau");
  auto g_d   = to_device<dview_1d>(g,   "g");
  auto w0_d  = to_device<dview_1d>(w0,  "w0");
  dview_1d ref_d("ref_d", n), tra_d("tra_d", n), ref2_d("ref2_d", n);
  Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,n), KOKKOS_LAMBDA (const int i) {
    F::two_stream(tau_d(i), g_d(i), w0_d(i), ref_d(i), tra_d(i));
    ref2_d(i) = F::two_stream_reflectance(tau_d(i), g_d(i), w0_d(i));
  });
  REQUIRE(count_diffs(ref_d,  ref_f)  == 0);
  REQUIRE(count_diffs(tra_d,  tra_f)  == 0);
  REQUIRE(count_diffs(ref2_d, ref2_f) == 0);
}

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_modis_bfb", "[cosp]")
{
  using C = F::C;
  constexpr int ntau = C::numMODISTauBins, npres = C::numMODISPresBins;

  for (const int ncol : {1, 3, 10}) {
    const int npoints = 300, nlev = 72;

    lview_2d phalf("phalf", npoints, nlev+1), pfull("pfull", npoints, nlev),
             at("at", npoints, nlev), qv("qv", npoints, nlev);
    lview_1d skt("skt", npoints);
    lview_3d tau("tau", npoints, ncol, nlev);
    make_thermo(phalf, pfull, at, qv, skt);
    make_dtau(tau);
    // Low clouds only in some points, so that the ISCCP cloud-top pressure gets used
    for (int j = 6; j < 60; ++j)
      for (int i = 0; i < ncol; ++i)
        for (int k = 0; k < nlev-12; ++k)
          tau(j,i,k) = 0;

    // Optical properties from the (already BFB) MODIS optics, with a mix of liquid, ice
    // and mixed-phase layers
    lview_3d tau_liq("tau_liq", npoints, ncol, nlev), tau_ice("tau_ice", npoints, ncol, nlev),
             size_liq("size_liq", npoints, ncol, nlev), size_ice("size_ice", npoints, ncol, nlev),
             liqfrac("liqfrac", npoints, ncol, nlev), g("g", npoints, ncol, nlev),
             w0("w0", npoints, ncol, nlev);
    for (int j = 0; j < npoints; ++j) {
      // Bias each point toward one phase so that all retrieved phases show up
      const double r = uniform(0,1);
      const double fmin = r < 0.4 ? 0.8 : (r < 0.7 ? 0.0 : 0.3);
      const double fmax = r < 0.4 ? 1.0 : (r < 0.7 ? 0.2 : 0.7);
      for (int i = 0; i < ncol; ++i) {
        for (int k = 0; k < nlev; ++k) {
          const double f = uniform(fmin, fmax);
          tau_liq(j,i,k)  = f*tau(j,i,k);
          tau_ice(j,i,k)  = tau(j,i,k) - tau_liq(j,i,k);
          size_liq(j,i,k) = uniform(2,35);
          size_ice(j,i,k) = uniform(3,100);
        }
      }
    }
    cosp_modis_optics_f(npoints, ncol, nlev, tau_liq.data(), size_liq.data(), tau_ice.data(),
                        size_ice.data(), liqfrac.data(), g.data(), w0.data());

    lview_2d isccp_ctp("isccp_ctp", npoints, ncol);
    for (int j = 0; j < npoints; ++j)
      for (int i = 0; i < ncol; ++i)
        isccp_ctp(j,i) = uniform(0,1) < 0.2 ? C::R_UNDEF : uniform(5000,100000);

    const double w0_cons = static_cast<double>(0.9999999f);
    // Edge cases: total optical depth right at the detection threshold, single scattering
    // albedo right at the (single precision) conservative scattering limit, and a liquid
    // fraction right at the phase threshold
    const double edge_tau[] = {0.3, std::nextafter(0.3, 0.0), 2.0, 2.0, 2.0, 2.0};
    const double edge_w0[]  = {-1, -1, w0_cons, std::nextafter(w0_cons, 1.0), 1.0, -1};
    for (int j = 0; j < 6; ++j) {
      for (int i = 0; i < ncol; ++i) {
        for (int k = 0; k < nlev; ++k) {
          tau(j,i,k) = k == nlev/2 || (j >= 2 && k == nlev/2+1) ? edge_tau[j] : 0.0;
          if (edge_w0[j] > 0) w0(j,i,k) = edge_w0[j];
          if (j == 5) liqfrac(j,i,k) = 0.7;
        }
      }
    }

    // Fortran: subcolumn, then column using the subcolumn outputs
    liview_2d phase_f("phase_f", npoints, ncol);
    lview_2d ctp_f("ctp_f", npoints, ncol), tau_ret_f("tau_ret_f", npoints, ncol),
             size_f("size_f", npoints, ncol);
    lview_3d fq_f("fq_f", npoints, ntau, npres);
    cosp_modis_subcolumn_f(npoints, ncol, nlev, phalf.data(), tau.data(), liqfrac.data(),
                           g.data(), w0.data(), isccp_ctp.data(), phase_f.data(), ctp_f.data(),
                           tau_ret_f.data(), size_f.data());
    cosp_modis_column_f(npoints, ncol, phase_f.data(), ctp_f.data(), tau_ret_f.data(),
                        size_f.data(), fq_f.data());

    // C++
    auto phalf_d   = to_device<dview_2d>(phalf,     "phalf");
    auto tau_d     = to_device<dview_3d>(tau,       "tau");
    auto liqfrac_d = to_device<dview_3d>(liqfrac,   "liqfrac");
    auto g_d       = to_device<dview_3d>(g,         "g");
    auto w0_d      = to_device<dview_3d>(w0,        "w0");
    auto ictp_d    = to_device<dview_2d>(isccp_ctp, "isccp_ctp");
    diview_2d phase_d("phase_d", npoints, ncol);
    dview_2d ctp_d("ctp_d", npoints, ncol), tau_ret_d("tau_ret_d", npoints, ncol),
             size_d("size_d", npoints, ncol);
    dview_3d fq_d("fq_d", npoints, ntau, npres);

    Kokkos::parallel_for(Policy2({0,0},{npoints,ncol}), KOKKOS_LAMBDA (const int j, const int i) {
      auto ph = Kokkos::subview(phalf_d,   j, Kokkos::ALL);
      auto t  = Kokkos::subview(tau_d,     j, i, Kokkos::ALL);
      auto lf = Kokkos::subview(liqfrac_d, j, i, Kokkos::ALL);
      auto gg = Kokkos::subview(g_d,       j, i, Kokkos::ALL);
      auto ww = Kokkos::subview(w0_d,      j, i, Kokkos::ALL);
      F::modis_subcolumn(nlev, ph, t, lf, gg, ww, ictp_d(j,i), phase_d(j,i), ctp_d(j,i),
                         tau_ret_d(j,i), size_d(j,i));
    });
    Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,npoints), KOKKOS_LAMBDA (const int j) {
      auto p  = Kokkos::subview(phase_d,   j, Kokkos::ALL);
      auto c  = Kokkos::subview(ctp_d,     j, Kokkos::ALL);
      auto t  = Kokkos::subview(tau_ret_d, j, Kokkos::ALL);
      auto s  = Kokkos::subview(size_d,    j, Kokkos::ALL);
      auto fq = Kokkos::subview(fq_d,      j, Kokkos::ALL, Kokkos::ALL);
      F::modis_column(ncol, p, c, t, s, fq);
    });

    REQUIRE(count_diffs(phase_d,   phase_f)   == 0);
    REQUIRE(count_diffs(ctp_d,     ctp_f)     == 0);
    REQUIRE(count_diffs(tau_ret_d, tau_ret_f) == 0);
    REQUIRE(count_diffs(size_d,    size_f)    == 0);
    REQUIRE(count_diffs(fq_d,      fq_f)      == 0);

    // Make sure the inputs exercise the interesting branches
    int nphase[4] = {0,0,0,0}, n_fill = 0, n_size = 0, n_isccp = 0, n_co2 = 0;
    double fq_total = 0;
    for (int j = 0; j < npoints; ++j) {
      for (int i = 0; i < ncol; ++i) {
        ++nphase[phase_f(j,i)];
        if (phase_f(j,i) == C::phaseIsNone) continue;
        if (size_f(j,i) == 1e-6*C::re_fill) ++n_fill;
        else ++n_size;
        if (ctp_f(j,i) == isccp_ctp(j,i)) ++n_isccp;
        else ++n_co2;
      }
      for (int a = 0; a < ntau; ++a)
        for (int b = 0; b < npres; ++b)
          fq_total += fq_f(j,a,b);
    }
    INFO("ncol=" << ncol << " none/liq/ice/undet=" << nphase[0] << "/" << nphase[1] << "/"
         << nphase[2] << "/" << nphase[3] << " n_fill=" << n_fill << " n_size=" << n_size
         << " n_isccp=" << n_isccp << " n_co2=" << n_co2);
    CHECK(nphase[0] > 0);
    CHECK(nphase[1] > 0);
    CHECK(nphase[2] > 0);
    CHECK(nphase[3] > 0);
    CHECK(n_fill > 0);
    CHECK(n_size > 0);
    CHECK(n_isccp > 0);
    CHECK(n_co2 > 0);
    CHECK(fq_total > 0);
    // The threshold edge cases: tau = 0.3 is detected, just below it is not
    CHECK(phase_f(0,0) != C::phaseIsNone);
    CHECK(phase_f(1,0) == C::phaseIsNone);
    CHECK(phase_f(5,0) == C::phaseIsLiquid);
  }
}
