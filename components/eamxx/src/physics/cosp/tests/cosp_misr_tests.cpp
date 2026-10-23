#include <catch2/catch.hpp>

#include "physics/cosp/tests/cosp_test_utils.hpp"

/*
 * BFB tests of the C++ MISR simulator, and of the COSP constants and bin tables, against
 * the Fortran COSP library.
 */

using namespace scream;
using namespace scream::cosp_test;

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_constants_bfb", "[cosp]")
{
  using C = F::C;
  double r_undef, misr_tauchk;
  lview_1d tau_bnds("tau_bnds", C::ntau+1), pres_bnds("pres_bnds", C::npres+1),
           hgt_bnds("hgt_bnds", C::nhgt+1);
  cosp_constants_f(&r_undef, tau_bnds.data(), pres_bnds.data(), hgt_bnds.data(), &misr_tauchk);

  REQUIRE(C::R_UNDEF == r_undef);
  REQUIRE(C::misr_tauchk == misr_tauchk);
  for (int i = 0; i <= C::ntau;  ++i) REQUIRE(C::tau_binBounds(i)  == tau_bnds(i));
  for (int i = 0; i <= C::npres; ++i) REQUIRE(C::pres_binBounds(i) == pres_bnds(i));
  for (int i = 0; i <= C::nhgt;  ++i) REQUIRE(C::hgt_binBounds(i)  == hgt_bnds(i));
}

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_misr_bfb", "[cosp]")
{
  using C = F::C;
  constexpr int ntau = C::numMISRTauBins, nhgt = C::numMISRHgtBins;

  for (const int ncol : {1, 3, 10}) {
    const int npoints = 200, nlev = 72;

    lview_3d dtau("dtau", npoints, ncol, nlev);
    lview_2d zfull("zfull", npoints, nlev), at("at", npoints, nlev);
    liview_1d sunlit("sunlit", npoints);
    make_profiles(zfull, at);
    make_dtau(dtau);
    for (int j = 0; j < npoints; ++j) sunlit(j) = uniform(0,1) < 0.75 ? 1 : 0;

    // Edge cases: single thin layers whose optical depth sits right at the thresholds,
    // which are single precision in the Fortran (0.2 vs 0.2f matters here)
    const double edge_taus[] = {0.2, 0.2000000001, 0.20000000298023224, 0.5, 0.5000000001,
                                C::misr_tauchk, std::nextafter(C::misr_tauchk, 1.0), 1.0e-8};
    int jedge = 0;
    for (const double tau_edge : edge_taus) {
      sunlit(jedge) = 1;
      for (int i = 0; i < ncol; ++i)
        for (int k = 0; k < nlev; ++k)
          dtau(jedge,i,k) = k == nlev/2 ? tau_edge : 0.0;
      ++jedge;
    }

    // Fortran: subcolumn, then column using the subcolumn outputs
    lview_2d tau_f("tau_f", npoints, ncol), ztop_f("ztop_f", npoints, ncol);
    lview_3d fq_f("fq_f", npoints, ntau, nhgt);
    cosp_misr_subcolumn_f(npoints, ncol, nlev, dtau.data(), zfull.data(), at.data(),
                          sunlit.data(), tau_f.data(), ztop_f.data());
    cosp_misr_column_f(npoints, ncol, ztop_f.data(), sunlit.data(), tau_f.data(), fq_f.data());

    // C++
    auto dtau_d  = to_device<dview_3d>(dtau,  "dtau");
    auto zfull_d = to_device<dview_2d>(zfull, "zfull");
    auto at_d    = to_device<dview_2d>(at,    "at");
    diview_1d sunlit_d("sunlit", npoints);
    {
      auto h = Kokkos::create_mirror_view(sunlit_d);
      Kokkos::deep_copy(h, sunlit);
      Kokkos::deep_copy(sunlit_d, h);
    }
    dview_2d tau_d("tau_d", npoints, ncol), ztop_d("ztop_d", npoints, ncol);
    dview_3d fq_d("fq_d", npoints, ntau, nhgt);

    Kokkos::parallel_for(Policy2({0,0},{npoints,ncol}), KOKKOS_LAMBDA (const int j, const int i) {
      auto dt = Kokkos::subview(dtau_d, j, i, Kokkos::ALL);
      auto z  = Kokkos::subview(zfull_d, j, Kokkos::ALL);
      auto t  = Kokkos::subview(at_d, j, Kokkos::ALL);
      F::misr_subcolumn(nlev, dt, z, t, sunlit_d(j), tau_d(j,i), ztop_d(j,i));
    });
    Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,npoints), KOKKOS_LAMBDA (const int j) {
      auto ztop = Kokkos::subview(ztop_d, j, Kokkos::ALL);
      auto tau  = Kokkos::subview(tau_d, j, Kokkos::ALL);
      auto fq   = Kokkos::subview(fq_d, j, Kokkos::ALL, Kokkos::ALL);
      F::misr_column(ncol, ztop, sunlit_d(j), tau, fq);
    });

    REQUIRE(count_diffs(tau_d,  tau_f)  == 0);
    REQUIRE(count_diffs(ztop_d, ztop_f) == 0);
    REQUIRE(count_diffs(fq_d,   fq_f)   == 0);

    // Make sure the inputs exercise the interesting branches: thin clouds flagged with
    // ztop = -1, clear subcolumns, interpolated tops, and a non-empty histogram
    int n_minus1 = 0, n_zero = 0, n_interp = 0;
    double fq_total = 0;
    for (int j = 0; j < npoints; ++j) {
      if (sunlit(j) != 1) continue;
      for (int i = 0; i < ncol; ++i) {
        if (ztop_f(j,i) == -1) ++n_minus1;
        else if (ztop_f(j,i) == 0) ++n_zero;
        else {
          bool on_level = false;
          for (int k = 0; k < nlev; ++k) on_level |= (ztop_f(j,i) == zfull(j,k));
          if (!on_level) ++n_interp;
        }
      }
      for (int a = 0; a < ntau; ++a)
        for (int b = 0; b < nhgt; ++b)
          fq_total += fq_f(j,a,b);
    }
    INFO("ncol=" << ncol << " n_minus1=" << n_minus1 << " n_zero=" << n_zero
         << " n_interp=" << n_interp);
    CHECK(n_minus1 > 0);
    CHECK(n_zero > 0);
    CHECK(n_interp > 0);
    CHECK(fq_total > 0);
  }
}
