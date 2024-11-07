#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestCloudSed : public UnitWrap::UnitTest<D>::Base {

void run_phys()
{
  // TODO
}

void run_bfb()
{
  auto engine = Base::get_engine();

  CloudSedData csds_baseline[] = {
    //         kts, kte, ktop, kbot, kdir,        dt,    inv_dt, do_predict_nc,     precip_liq_surf,
    CloudSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,         false,     0.0),
    CloudSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,         false,     0.0),
    CloudSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,          true,     0.0),
    CloudSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,          true,     0.0),
    CloudSedData(1,  72,   27,   27,   -1, 1.800E+03, 5.556E-04,          true,     0.0),
  };

  static constexpr Int num_runs = sizeof(csds_baseline) / sizeof(CloudSedData);

  // Set up random input data
  for (auto& d : csds_baseline) {
    d.randomize(engine, { {d.qc_incld, {C::QSMALL/2, C::QSMALL*2}} });
  }

  // Create copies of data for use by cxx. Needs to happen before reads so that
  // inout data is in original state
  CloudSedData csds_cxx[num_runs] = {
    CloudSedData(csds_baseline[0]),
    CloudSedData(csds_baseline[1]),
    CloudSedData(csds_baseline[2]),
    CloudSedData(csds_baseline[3]),
    CloudSedData(csds_baseline[4]),
  };

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (auto& d : csds_baseline) {
      d.read(Base::m_fid);
    }
  }

  // Get data from cxx
  for (auto& d : csds_cxx) {
    cloud_sedimentation_host(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                          d.qc_incld, d.rho, d.inv_rho, d.cld_frac_l, d.acn, d.inv_dz,
                          d.dt, d.inv_dt, d.do_predict_nc,
                          d.qc, d.nc, d.nc_incld, d.mu_c, d.lamc, &d.precip_liq_surf, d.qc_tend, d.nc_tend);
  }

  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < num_runs; ++i) {
      // Due to pack issues, we must restrict checks to the active k space
      Int start = std::min(csds_baseline[i].kbot, csds_baseline[i].ktop) - 1; // 0-based indx
      Int end   = std::max(csds_baseline[i].kbot, csds_baseline[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(csds_baseline[i].qc[k]       == csds_cxx[i].qc[k]);
        REQUIRE(csds_baseline[i].nc[k]       == csds_cxx[i].nc[k]);
        REQUIRE(csds_baseline[i].nc_incld[k] == csds_cxx[i].nc_incld[k]);
        REQUIRE(csds_baseline[i].mu_c[k]     == csds_cxx[i].mu_c[k]);
        REQUIRE(csds_baseline[i].lamc[k]     == csds_cxx[i].lamc[k]);
        REQUIRE(csds_baseline[i].qc_tend[k]  == csds_cxx[i].qc_tend[k]);
        REQUIRE(csds_baseline[i].nc_tend[k]  == csds_cxx[i].nc_tend[k]);
      }
      REQUIRE(csds_baseline[i].precip_liq_surf == csds_cxx[i].precip_liq_surf);
    }
  }
  else if (this->m_baseline_action == GENERATE) {
    for (Int i = 0; i < num_runs; ++i) {
      csds_cxx[i].write(Base::m_fid);
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
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCloudSed;

  T t;
  t.run_phys();
  t.run_bfb();
}

} // namespace
