#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

#include <iostream>
#include <cstdlib>
#include <unistd.h>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPblintdHeight : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
    static const auto approx_zero = Approx(0.0).margin(1e-16);
    static constexpr Int shcol = 4;
    static constexpr Int nlev = 9;

    // Tests for the subroutine pblintd_height
    // Perform a series of tests to ensure that subroutine returns values
    //  as expected

    // TEST ONE
    // Given input where the value of ustar is increasing, verify that
    //  the PBL depth also increases.  Input temperature should be stable

    // Define the midpoint heights [m]
    static constexpr Real z[nlev] = {5000, 4000, 3000, 2000, 1000, 750, 500, 250, 100};
    // Define the u wind [m/s]
    static constexpr Real u[nlev] = {10, 10, 10, 9, 9, 8, 7, 6, 5};
    // Define the v wind [m/s]
    static constexpr Real v[nlev] = {-10, -10, -10, -9, -9, -8, -7, -6, -5};
    // Define the virtual potential temperature [K]
    static constexpr Real thv[nlev] = {320, 315, 314, 312, 311, 310, 302, 302, 300};
    // Define the surface friction velocity [m4/s3]
    Real ustar = ustar_min;

    // Initialize dtata structure for bridging to F90
    PblintdHeightData SDS(shcol,nlev,nlev);

    // Require that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    // This test requires more than one column
    REQUIRE(shcol > 1);

    // Fill in test data on the zt_grid
    for(Int s = 0; s < shcol; ++s) {
      SDS.ustar[s] = ustar + s; // Add to ustar with each column
      SDS.check[s] = true;
      SDS.thv_ref[s] = thv[nlev-1];
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.z[offset] = z[n];
        SDS.u[offset] = u[n];
        SDS.v[offset] = v[n];
        SDS.thv[offset] = thv[n];
        // Should be initialized to zero everywhere
        SDS.rino[offset] = 0;
      }
    }

    // check to make sure the input data makes sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.ustar[s] >= ustar_min);
      REQUIRE(SDS.check[s] == true);
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.z[offset] > 0);
        REQUIRE(std::abs(SDS.u[offset]) < 100);
        REQUIRE(std::abs(SDS.v[offset]) < 100);
        REQUIRE(SDS.thv[offset] < 1000);
        REQUIRE(SDS.thv[offset] > 150);
      }
    }

    // check that as column increases then does ustar
    for(Int s = 0; s < shcol-1; ++s) {
      REQUIRE(SDS.ustar[s+1] > SDS.ustar[s]);
    }

    // Call the C++ implementation
    pblintd_height(SDS);

    // Check the result
    for(Int s = 0; s < shcol; ++s) {
      // Verify PBL height falls within reasonable bounds
      REQUIRE(SDS.pblh[s] > z[nlev-1]); // should be higher than lowest grid
      REQUIRE(SDS.pblh[s] < z[0]); // and lower than highest grid
      REQUIRE(SDS.check[s] == false);
    }

    for(Int s = 0; s < shcol-1; ++s) {
      // Verify that as ustar increases, PBLH increases
      REQUIRE(SDS.pblh[s+1] > SDS.pblh[s]);
    }

    // Save result for checking with next test
    Real pblh_save[shcol];
    for(Int s = 0; s < shcol; ++s) {
      pblh_save[s] = SDS.pblh[s];
    }

    // TEST TWO
    // Now use the same inputs as the last test BUT modify the surface
    //  temperature value so that it is unstable.  The PBLH should be
    //  larger everywhere.

    // Replace lowest input temperature to a value that is unstable
    for(Int s = 0; s < shcol; ++s) {
      SDS.check[s] = true; // reinitialize
      const auto offset = nlev-1 + s * nlev;
      SDS.thv[offset] = SDS.thv[offset-1] + 2;
    }

    // Call the C++ implementation
    pblintd_height(SDS);

    // Check the result
    for(Int s = 0; s < shcol; ++s) {
      // Verify PBL height falls within reasonable bounds
      REQUIRE(SDS.pblh[s] > z[nlev-1]); // should be higher than lowest grid
      REQUIRE(SDS.pblh[s] < z[0]); // and lower than highest grid
      REQUIRE(SDS.check[s] == false);
    }

    // Verify that the PBLh with unstable atmosphere is LARGER everywhere
    //  than the tests with conditionally stable atmosphere
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > pblh_save[s]);
    }

    // TEST THREE
    // Given profiles where all variables are homogenous, verify that
    //  PBLH everywhere is returned as zero (in real code the subroutine
    //  will be called again to ammeliorate this, tested in upper level function
    //  for the PBLH).

    // Fill in test data on the zt_grid, all other data recycled
    for(Int s = 0; s < shcol; ++s) {
      SDS.check[s] = true;
      SDS.thv_ref[s] = 300;
      SDS.pblh[s] = 0;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Make all input, homogenous
        SDS.u[offset] = 5;
        SDS.v[offset] = -2;
        SDS.thv[offset] = 300;
        // Should be initialized to zero everywhere
        SDS.rino[offset] = 0;
      }
    }

    // Call the C++ implementation
    pblintd_height(SDS);

    // Check that PBLH is zero (not modified) everywhere
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] == approx_zero);
    }

  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // [HUNT] print the raw libc rand() value, to see if it varies from run to run
    const int hunt_rand_raw = rand();
    Int npbl_rand = hunt_rand_raw%72 + 1;
    {
      char hn[256] = "?"; gethostname(hn, sizeof(hn));
      const char* omp = getenv("OMP_NUM_THREADS");
      std::cout << "[HUNT] pblintd_height_bfb: mode=" << (this->m_baseline_action==COMPARE ? "COMPARE" : (this->m_baseline_action==GENERATE ? "GENERATE" : "NONE"))
                << " rand()=" << hunt_rand_raw << " RAND_MAX=" << RAND_MAX
                << " npbl_rand=" << npbl_rand << " host=" << hn
                << " OMP_NUM_THREADS=" << (omp ? omp : "unset") << std::endl;
    }

    PblintdHeightData baseline_data[] = {
      PblintdHeightData(10, 72, 1),
      PblintdHeightData(10, 72, 72),
      PblintdHeightData(10, 72, npbl_rand),
      PblintdHeightData(10, 12, 1),
      PblintdHeightData(7, 16, 1),
      PblintdHeightData(2, 7, 1),
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    PblintdHeightData cxx_data[] = {
      PblintdHeightData(baseline_data[0]),
      PblintdHeightData(baseline_data[1]),
      PblintdHeightData(baseline_data[2]),
      PblintdHeightData(baseline_data[3]),
      PblintdHeightData(baseline_data[4]),
      PblintdHeightData(baseline_data[5]),
    };

    // [HUNT] pristine copies of the (regenerated) inputs, to rerun the kernel on identical inputs
    PblintdHeightData orig_data[] = {
      PblintdHeightData(baseline_data[0]),
      PblintdHeightData(baseline_data[1]),
      PblintdHeightData(baseline_data[2]),
      PblintdHeightData(baseline_data[3]),
      PblintdHeightData(baseline_data[4]),
      PblintdHeightData(baseline_data[5]),
    };

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      pblintd_height(d);
    }

    // [HUNT] diagnostics (host side only). Note: the baseline file stores the data ARRAYS only, not
    // the scalar dims (shcol, nlev, npbl). So after read(), baseline_data holds the arrays from the
    // generation run (inputs and outputs), while cxx_data/orig_data hold the regenerated inputs.
    if (this->m_baseline_action == COMPARE) {
      auto ndiff = [](const Real* a, const Real* b, Int n) {
        Int nd = 0;
        for (Int k=0; k<n; ++k) { if (a[k]!=b[k]) ++nd; }
        return nd;
      };
      auto ndiff_b = [](const bool* a, const bool* b, Int n) {
        Int nd = 0;
        for (Int k=0; k<n; ++k) { if (a[k]!=b[k]) ++nd; }
        return nd;
      };
      static constexpr Int hunt_nsets = sizeof(baseline_data) / sizeof(PblintdHeightData);
      for (Int i = 0; i < hunt_nsets; ++i) {
        auto& b = baseline_data[i];
        auto& c = cxx_data[i];
        const Int nc = c.total(c.pblh);
        const Int nl = c.total(c.z);
        // Inputs: file vs regenerated. If these differ, the inputs generation (engine/seed) differs.
        const Int d_in = ndiff(b.z,c.z,nl) + ndiff(b.u,c.u,nl) + ndiff(b.v,c.v,nl) + ndiff(b.thv,c.thv,nl)
                       + ndiff(b.ustar,c.ustar,nc) + ndiff(b.thv_ref,c.thv_ref,nc);
        // Outputs: file vs cxx
        const Int d_pblh  = ndiff(b.pblh,c.pblh,nc);
        const Int d_check = ndiff_b(b.check,c.check,nc);
        const Int d_rino  = ndiff(b.rino,c.rino,nl);
        // Determinism: rerun the kernel (same npbl) on identical inputs a few times
        Int d_rerun = 0;
        for (int r=0; r<5; ++r) {
          PblintdHeightData t(orig_data[i]);
          pblintd_height(t);
          d_rerun += ndiff(t.pblh,c.pblh,nc) + ndiff_b(t.check,c.check,nc) + ndiff(t.rino,c.rino,nl);
        }
        std::cout << "[HUNT] set " << i << " (shcol,nlev,npbl_cxx)=(" << c.shcol << "," << c.nlev << "," << c.npbl << ")"
                  << " input_diffs(file vs regen)=" << d_in
                  << " pblh_diffs=" << d_pblh << " check_diffs=" << d_check << " rino_diffs=" << d_rino
                  << " kernel_rerun_diffs(5 reruns)=" << d_rerun << std::endl;
        if (d_pblh>0 || d_check>0) {
          // Which npbl reproduces the baseline outputs, given the regenerated inputs?
          std::cout << "[HUNT]   set " << i << ": npbl values reproducing the baseline pblh/check:";
          int nmatch = 0;
          for (Int np=1; np<=c.nlev; ++np) {
            PblintdHeightData t(orig_data[i]);
            t.npbl = np;
            pblintd_height(t);
            if (ndiff(t.pblh,b.pblh,nc)==0 && ndiff_b(t.check,b.check,nc)==0) {
              std::cout << " " << np;
              ++nmatch;
            }
          }
          std::cout << (nmatch==0 ? " NONE" : "") << std::endl;
          for (Int k=0; k<nc && k<3; ++k) {
            std::cout << "[HUNT]   set " << i << " col " << k << ": ustar=" << c.ustar[k] << " thv_ref=" << c.thv_ref[k]
                      << " pblh(file)=" << b.pblh[k] << " pblh(cxx)=" << c.pblh[k]
                      << " check(file)=" << b.check[k] << " check(cxx)=" << c.check[k] << std::endl;
          }
        }
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(PblintdHeightData);
      for (Int i = 0; i < num_runs; ++i) {
        PblintdHeightData& d_baseline = baseline_data[i];
        PblintdHeightData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.pblh); ++k) {
          REQUIRE(d_baseline.total(d_baseline.pblh) == d_cxx.total(d_cxx.pblh));
          REQUIRE(d_baseline.pblh[k] == d_cxx.pblh[k]);
          REQUIRE(d_baseline.total(d_baseline.pblh) == d_cxx.total(d_cxx.check));
          REQUIRE(d_baseline.check[k] == d_cxx.check[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_ofile);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("pblintd_height_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdHeight;

  TestStruct().run_property();
}

TEST_CASE("pblintd_height_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdHeight;

  TestStruct().run_bfb();
}

} // empty namespace
