#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

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
  const std::array< std::pair<Real, Real>, CloudSedData::NUM_ARRAYS > ranges = {
    std::make_pair(5.100E-03, 9.952E-07), // qc_incld_range
    std::make_pair(4.056E-03, 1.153E+00), // rho_range
    std::make_pair(0,         1.),        // inv_rho (ignored)
    std::make_pair(0.9, 1.1),             // cld_frac_l_range
    std::make_pair(2.959E+07, 5.348E+07), // acn_range
    std::make_pair(2.863E-05, 8.141E-03), // inv_dz_range
    std::make_pair(7.701E-16, 2.119E-04), // qc_range
    std::make_pair(7.701E-16, 2.119E-04), // nc_range
    std::make_pair(9.952E+05, 1.734E+08), // nc_incld_range
    std::make_pair(5.722E+00, 1.253E+01), // mu_c_range
    std::make_pair(3.381E+05, 2.519E+06), // lamc_range
    std::make_pair(9.952E-07, 9.982E-07), // qc_tend_range
    std::make_pair(9.952E+05, 1.743E+08), // nc_tend_range
  };

  CloudSedData csds_fortran[] = {
    //         kts, kte, ktop, kbot, kdir,        dt,       inv_dt, do_predict_nc, precip_liq_surf, ranges
    CloudSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,         false,     0.0, ranges),
    CloudSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,         false,     0.0, ranges),
    CloudSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,          true,     0.0, ranges),
    CloudSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,          true,     0.0, ranges),
    CloudSedData(1,  72,   27,   27,   -1, 1.800E+03, 5.556E-04,          true,     0.0, ranges),
  };

  static constexpr Int num_runs = sizeof(csds_fortran) / sizeof(CloudSedData);

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
  for (Int i = 0; i < num_runs; ++i) {
    cloud_sedimentation(csds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    CloudSedData& d = csds_cxx[i];
    cloud_sedimentation_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                          d.qc_incld, d.rho, d.inv_rho, d.cld_frac_l, d.acn, d.inv_dz,
                          d.dt, d.inv_dt, d.do_predict_nc,
                          d.qc, d.nc, d.nc_incld, d.mu_c, d.lamc, &d.precip_liq_surf, d.qc_tend, d.nc_tend);
  }

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
