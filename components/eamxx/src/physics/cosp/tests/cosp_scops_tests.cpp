#include <catch2/catch.hpp>

#include "physics/cosp/tests/cosp_test_utils.hpp"

/*
 * BFB tests of the C++ subcolumn generator (SCOPS) against the Fortran COSP library.
 */

using namespace scream;
using namespace scream::cosp_test;

// -------------------------------------------------------------------------------------------
TEST_CASE("cosp_scops_bfb", "[cosp]")
{
  using C = F::C;
  for (const int overlap : {1, 2, 3, 7}) {
    for (const int ncol : {2, 10, 50}) {
      const int npoints = 200, nlev = 72;
      lview_2d cc("cc", npoints, nlev), conv("conv", npoints, nlev);
      liview_1d seeds("seeds", npoints);
      for (int j = 0; j < npoints; ++j) {
        // Seeds as subsample_and_optics makes them, plus a few large ones
        seeds(j) = j < 5 ? int(uniform(-2.0e9, 2.0e9)) : int(uniform(0,1)*100000) + 1;
        const double cld_prob = uniform(0,1);
        const bool with_conv = uniform(0,1) < 0.3;
        for (int k = 0; k < nlev; ++k) {
          const double r = uniform(0,1);
          cc(j,k) = r < cld_prob ? (r < 0.2*cld_prob ? 1.0 : uniform(0,1)) : 0.0;
          conv(j,k) = with_conv && uniform(0,1) < 0.3 ? uniform(0,1)*cc(j,k) : 0.0;
          // Edge cases: cloud fractions equal to the subcolumn positions in double
          // precision, while the Fortran computes those in single precision
          if (j >= 5 && j < 10) {
            cc(j,k) = (k % ncol + 0.5)/ncol;
            if (j >= 7) conv(j,k) = cc(j,k);
          }
        }
      }

      lview_3d fo_f("fo_f", npoints, ncol, nlev);
      cosp_scops_f(npoints, nlev, ncol, seeds.data(), cc.data(), conv.data(), overlap,
                   fo_f.data());

      auto cc_d   = to_device<dview_2d>(cc,   "cc");
      auto conv_d = to_device<dview_2d>(conv, "conv");
      diview_1d seeds_d("seeds", npoints);
      {
        auto h = Kokkos::create_mirror_view(seeds_d);
        Kokkos::deep_copy(h, seeds);
        Kokkos::deep_copy(seeds_d, h);
      }
      dview_3d fo_d("fo_d", npoints, ncol, nlev);
      dview_2d thr_d("thr_d", npoints, ncol);
      Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,npoints), KOKKOS_LAMBDA (const int j) {
        F::scops(nlev, ncol, seeds_d(j), seeds_d(npoints-1),
                 Kokkos::subview(cc_d, j, Kokkos::ALL), Kokkos::subview(conv_d, j, Kokkos::ALL),
                 overlap, Kokkos::subview(fo_d, j, Kokkos::ALL, Kokkos::ALL),
                 Kokkos::subview(thr_d, j, Kokkos::ALL));
      });

      INFO("overlap=" << overlap << " ncol=" << ncol);
      REQUIRE(count_diffs(fo_d, fo_f) == 0);

      // All three subcolumn types show up
      int n[3] = {0,0,0};
      for (int j = 0; j < npoints; ++j)
        for (int i = 0; i < ncol; ++i)
          for (int k = 0; k < nlev; ++k)
            ++n[int(fo_f(j,i,k))];
      CHECK(n[0] > 0);
      CHECK(n[C::I_LSC] > 0);
      CHECK(n[C::I_CVC] > 0);
    }
  }
}
