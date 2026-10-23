#include <catch2/catch.hpp>

#include "physics/cosp/tests/cosp_test_utils.hpp"

/*
 * BFB tests of the C++ ISCCP simulator (ICARUS) against the Fortran COSP library.
 */

using namespace scream;
using namespace scream::cosp_test;

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_isccp_bfb", "[cosp]")
{
  using C = F::C;
  constexpr int ntau = C::numISCCPTauBins, npres = C::numISCCPPresBins;

  double tauchk_f;
  cosp_isccp_constants_f(&tauchk_f);
  REQUIRE(C::isccp_tauchk == tauchk_f);

  const double emsfc_lw = 0.99;
  for (const int ncol : {1, 3, 10}) {
    const int npoints = 300, nlev = 72;

    lview_2d phalf("phalf", npoints, nlev+1), pfull("pfull", npoints, nlev),
             at("at", npoints, nlev), qv("qv", npoints, nlev);
    lview_1d skt("skt", npoints);
    lview_3d dtau("dtau", npoints, ncol, nlev), dem("dem", npoints, ncol, nlev);
    liview_1d sunlit("sunlit", npoints);
    make_thermo(phalf, pfull, at, qv, skt);
    make_dtau(dtau);
    for (int j = 0; j < npoints; ++j) {
      sunlit(j) = uniform(0,1) < 0.8 ? 1 : 0;
      for (int i = 0; i < ncol; ++i)
        for (int k = 0; k < nlev; ++k)
          dem(j,i,k) = dtau(j,i,k) > 0 ? 1 - std::exp(-0.5*dtau(j,i,k)) : 0.0;
    }

    // Fortran
    lview_2d boxtau_f("boxtau_f", npoints, ncol), boxptop_f("boxptop_f", npoints, ncol),
             boxttop_f("boxttop_f", npoints, ncol);
    lview_3d fq_f("fq_f", npoints, ntau, npres);
    lview_1d cldarea_f("cldarea_f", npoints);
    cosp_icarus_f(npoints, ncol, nlev, sunlit.data(), dtau.data(), dem.data(), skt.data(),
                  emsfc_lw, qv.data(), at.data(), pfull.data(), phalf.data(), boxtau_f.data(),
                  boxptop_f.data(), boxttop_f.data(), fq_f.data(), cldarea_f.data());

    // C++
    auto phalf_d = to_device<dview_2d>(phalf, "phalf");
    auto pfull_d = to_device<dview_2d>(pfull, "pfull");
    auto at_d    = to_device<dview_2d>(at,    "at");
    auto qv_d    = to_device<dview_2d>(qv,    "qv");
    auto skt_d   = to_device<dview_1d>(skt,   "skt");
    auto dtau_d  = to_device<dview_3d>(dtau,  "dtau");
    auto dem_d   = to_device<dview_3d>(dem,   "dem");
    diview_1d sunlit_d("sunlit", npoints);
    {
      auto h = Kokkos::create_mirror_view(sunlit_d);
      Kokkos::deep_copy(h, sunlit);
      Kokkos::deep_copy(sunlit_d, h);
    }
    dview_2d boxtau_d("boxtau_d", npoints, ncol), boxptop_d("boxptop_d", npoints, ncol),
             boxttop_d("boxttop_d", npoints, ncol);
    dview_3d fq_d("fq_d", npoints, ntau, npres);
    dview_1d cldarea_d("cldarea_d", npoints);

    Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,npoints), KOKKOS_LAMBDA (const int j) {
      auto q  = Kokkos::subview(qv_d,    j, Kokkos::ALL);
      auto t  = Kokkos::subview(at_d,    j, Kokkos::ALL);
      auto pf = Kokkos::subview(pfull_d, j, Kokkos::ALL);
      auto ph = Kokkos::subview(phalf_d, j, Kokkos::ALL);
      const auto pt = F::icarus_point(nlev, skt_d(j), emsfc_lw, q, t, pf, ph);
      for (int i = 0; i < ncol; ++i) {
        auto dt = Kokkos::subview(dtau_d, j, i, Kokkos::ALL);
        auto de = Kokkos::subview(dem_d,  j, i, Kokkos::ALL);
        F::icarus_subcolumn(nlev, sunlit_d(j), pt, dt, de, skt_d(j), emsfc_lw, q, t, pf, ph,
                            boxtau_d(j,i), boxptop_d(j,i), boxttop_d(j,i));
      }
      auto bt = Kokkos::subview(boxtau_d,  j, Kokkos::ALL);
      auto bp = Kokkos::subview(boxptop_d, j, Kokkos::ALL);
      auto fq = Kokkos::subview(fq_d, j, Kokkos::ALL, Kokkos::ALL);
      F::icarus_column(ncol, sunlit_d(j), bt, bp, fq, cldarea_d(j));
    });

    COSP_REQUIRE_MATCH(boxtau_d, boxtau_f);
    COSP_REQUIRE_MATCH(boxptop_d, boxptop_f);
    COSP_REQUIRE_MATCH(boxttop_d, boxttop_f);
    COSP_REQUIRE_MATCH(fq_d, fq_f);
    COSP_REQUIRE_MATCH(cldarea_d, cldarea_f);

    // Coverage: cloudy subcolumns with valid tops, at a range of pressures
    int ncloudy = 0, nhigh = 0, nlow = 0;
    for (int j = 0; j < npoints; ++j) {
      for (int i = 0; i < ncol; ++i) {
        if (boxtau_f(j,i) > 0) {
          ++ncloudy;
          if (boxptop_f(j,i) < 44000) ++nhigh;
          if (boxptop_f(j,i) > 68000) ++nlow;
        }
      }
    }
    INFO("ncol=" << ncol << " ncloudy=" << ncloudy << " nhigh=" << nhigh << " nlow=" << nlow);
    CHECK(ncloudy > npoints*ncol/10);
    CHECK(nhigh > 0);
    CHECK(nlow > 0);
  }
}
