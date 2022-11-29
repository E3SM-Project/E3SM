#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestCloudSed {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  auto engine = setup_random_test();

  CloudSedData csds_fortran[] = {
    //         kts, kte, ktop, kbot, kdir,        dt,    inv_dt, do_predict_nc,     precip_liq_surf,
    CloudSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,         false,     0.0),
    CloudSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,         false,     0.0),
    CloudSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,          true,     0.0),
    CloudSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,          true,     0.0),
    CloudSedData(1,  72,   27,   27,   -1, 1.800E+03, 5.556E-04,          true,     0.0),
  };

  static constexpr Int num_runs = sizeof(csds_fortran) / sizeof(CloudSedData);

  // Set up random input data
  for (auto& d : csds_fortran) {
    d.randomize(engine, { {d.qc_incld, {C::QSMALL/2, C::QSMALL*2}} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  CloudSedData csds_cxx[num_runs] = {
    CloudSedData(csds_fortran[0]),
    CloudSedData(csds_fortran[1]),
    CloudSedData(csds_fortran[2]),
    CloudSedData(csds_fortran[3]),
    CloudSedData(csds_fortran[4]),
  };

  // Get data from fortran
  for (auto& d : csds_fortran) {
    cloud_sedimentation(d);
  }

  // Get data from cxx
  for (auto& d : csds_cxx) {
    cloud_sedimentation_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                          d.qc_incld, d.rho, d.inv_rho, d.cld_frac_l, d.acn, d.inv_dz,
                          d.dt, d.inv_dt, d.do_predict_nc,
                          d.qc, d.nc, d.nc_incld, d.mu_c, d.lamc, &d.precip_liq_surf, d.qc_tend, d.nc_tend);
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      // Due to pack issues, we must restrict checks to the active k space
      Int start = std::min(csds_fortran[i].kbot, csds_fortran[i].ktop) - 1; // 0-based indx
      Int end   = std::max(csds_fortran[i].kbot, csds_fortran[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(csds_fortran[i].qc[k]       == csds_cxx[i].qc[k]);
        REQUIRE(csds_fortran[i].nc[k]       == csds_cxx[i].nc[k]);
        REQUIRE(csds_fortran[i].nc_incld[k] == csds_cxx[i].nc_incld[k]);
        REQUIRE(csds_fortran[i].mu_c[k]     == csds_cxx[i].mu_c[k]);
        REQUIRE(csds_fortran[i].lamc[k]     == csds_cxx[i].lamc[k]);
        REQUIRE(csds_fortran[i].qc_tend[k]  == csds_cxx[i].qc_tend[k]);
        REQUIRE(csds_fortran[i].nc_tend[k]  == csds_cxx[i].nc_tend[k]);
      }
      REQUIRE(csds_fortran[i].precip_liq_surf == csds_cxx[i].precip_liq_surf);
    }
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_cloud_sed", "[p3_functions]")
{
  using TCS = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCloudSed;

  TCS::run_phys();
  TCS::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
