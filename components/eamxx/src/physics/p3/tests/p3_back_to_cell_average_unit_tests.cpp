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
struct UnitWrap::UnitTest<D>::TestBackToCellAverage {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  auto engine = setup_random_test();

  // Generate n test structs, each populated with random data (values within
  // [0,1]) by the default constructor.
  BackToCellAverageData back_to_cell_average_data[Spack::n];
  for (Int i = 0; i < Spack::n; ++i) {
    back_to_cell_average_data[i].randomize(engine);
  }

  // Sync to device.
  view_1d<BackToCellAverageData> device_data("back_to_cell_average", Spack::n);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&back_to_cell_average_data[0], &back_to_cell_average_data[0] + Spack::n,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < Spack::n; ++i) {
    back_to_cell_average(back_to_cell_average_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack cld_frac_l, cld_frac_r, cld_frac_i, qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend, nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend,
      nr_evap_tend, ncautr, qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend, qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend,
      qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend,
      nr_collect_tend, ni_selfcollect_tend, qv2qi_vapdep_tend, nr2ni_immers_freeze_tend,
      ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend, qc2qi_berg_tend;
    for (Int s = 0; s < Spack::n; ++s) {
      cld_frac_l[s]               = device_data[s].cld_frac_l;
      cld_frac_r[s]               = device_data[s].cld_frac_r;
      cld_frac_i[s]               = device_data[s].cld_frac_i;
      qc2qr_accret_tend[s]        = device_data[s].qc2qr_accret_tend;
      qr2qv_evap_tend[s]          = device_data[s].qr2qv_evap_tend;
      qc2qr_autoconv_tend[s]      = device_data[s].qc2qr_autoconv_tend;
      nc_accret_tend[s]           = device_data[s].nc_accret_tend;
      nc_selfcollect_tend[s]      = device_data[s].nc_selfcollect_tend;
      nc2nr_autoconv_tend[s]      = device_data[s].nc2nr_autoconv_tend;
      nr_selfcollect_tend[s]      = device_data[s].nr_selfcollect_tend;
      nr_evap_tend[s]             = device_data[s].nr_evap_tend;
      ncautr[s]                   = device_data[s].ncautr;
      qi2qv_sublim_tend[s]        = device_data[s].qi2qv_sublim_tend;
      nr_ice_shed_tend[s]         = device_data[s].nr_ice_shed_tend;
      qc2qi_hetero_freeze_tend[s] = device_data[s].qc2qi_hetero_freeze_tend;
      qr2qi_collect_tend[s]       = device_data[s].qr2qi_collect_tend;
      qc2qr_ice_shed_tend[s]      = device_data[s].qc2qr_ice_shed_tend;
      qi2qr_melt_tend[s]          = device_data[s].qi2qr_melt_tend;
      qc2qi_collect_tend[s]       = device_data[s].qc2qi_collect_tend;
      qr2qi_immers_freeze_tend[s] = device_data[s].qr2qi_immers_freeze_tend;
      ni2nr_melt_tend[s]          = device_data[s].ni2nr_melt_tend;
      nc_collect_tend[s]          = device_data[s].nc_collect_tend;
      ncshdc[s]                   = device_data[s].ncshdc;
      nc2ni_immers_freeze_tend[s] = device_data[s].nc2ni_immers_freeze_tend;
      nr_collect_tend[s]          = device_data[s].nr_collect_tend;
      ni_selfcollect_tend[s]      = device_data[s].ni_selfcollect_tend;
      qv2qi_vapdep_tend[s]        = device_data[s].qv2qi_vapdep_tend;
      nr2ni_immers_freeze_tend[s] = device_data[s].nr2ni_immers_freeze_tend;
      ni_sublim_tend[s]           = device_data[s].ni_sublim_tend;
      qv2qi_nucleat_tend[s]       = device_data[s].qv2qi_nucleat_tend;
      ni_nucleat_tend[s]          = device_data[s].ni_nucleat_tend;
      qc2qi_berg_tend[s]          = device_data[s].qc2qi_berg_tend;
    }

    Functions::back_to_cell_average(cld_frac_l, cld_frac_r, cld_frac_i, qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend,
      nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr, qi2qv_sublim_tend, nr_ice_shed_tend,
      qc2qi_hetero_freeze_tend, qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend,
      qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend,
      nr_collect_tend, ni_selfcollect_tend, qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend, qc2qi_berg_tend);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      device_data(s).qc2qr_accret_tend        = qc2qr_accret_tend[s];
      device_data(s).qr2qv_evap_tend          = qr2qv_evap_tend[s];
      device_data(s).qc2qr_autoconv_tend      = qc2qr_autoconv_tend[s];
      device_data(s).nc_accret_tend           = nc_accret_tend[s];
      device_data(s).nc_selfcollect_tend      = nc_selfcollect_tend[s];
      device_data(s).nc2nr_autoconv_tend      = nc2nr_autoconv_tend[s];
      device_data(s).nr_selfcollect_tend      = nr_selfcollect_tend[s];
      device_data(s).nr_evap_tend             = nr_evap_tend[s];
      device_data(s).ncautr                   = ncautr[s];
      device_data(s).qi2qv_sublim_tend        = qi2qv_sublim_tend[s];
      device_data(s).nr_ice_shed_tend         = nr_ice_shed_tend[s];
      device_data(s).qc2qi_hetero_freeze_tend = qc2qi_hetero_freeze_tend[s];
      device_data(s).qr2qi_collect_tend       = qr2qi_collect_tend[s];
      device_data(s).qc2qr_ice_shed_tend      = qc2qr_ice_shed_tend[s];
      device_data(s).qi2qr_melt_tend          = qi2qr_melt_tend[s];
      device_data(s).qc2qi_collect_tend       = qc2qi_collect_tend[s];
      device_data(s).qr2qi_immers_freeze_tend = qr2qi_immers_freeze_tend[s];
      device_data(s).ni2nr_melt_tend          = ni2nr_melt_tend[s];
      device_data(s).nc_collect_tend          = nc_collect_tend[s];
      device_data(s).ncshdc                   = ncshdc[s];
      device_data(s).nc2ni_immers_freeze_tend = nc2ni_immers_freeze_tend[s];
      device_data(s).nr_collect_tend          = nr_collect_tend[s];
      device_data(s).ni_selfcollect_tend      = ni_selfcollect_tend[s];
      device_data(s).qv2qi_vapdep_tend        = qv2qi_vapdep_tend[s];
      device_data(s).nr2ni_immers_freeze_tend = nr2ni_immers_freeze_tend[s];
      device_data(s).ni_sublim_tend           = ni_sublim_tend[s];
      device_data(s).qv2qi_nucleat_tend       = qv2qi_nucleat_tend[s];
      device_data(s).ni_nucleat_tend          = ni_nucleat_tend[s];
      device_data(s).qc2qi_berg_tend          = qc2qi_berg_tend[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  if (SCREAM_BFB_TESTING) {
    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(back_to_cell_average_data[s].qc2qr_accret_tend        == host_data[s].qc2qr_accret_tend);
      REQUIRE(back_to_cell_average_data[s].qr2qv_evap_tend          == host_data[s].qr2qv_evap_tend);
      REQUIRE(back_to_cell_average_data[s].qc2qr_autoconv_tend      == host_data[s].qc2qr_autoconv_tend);
      REQUIRE(back_to_cell_average_data[s].nc_accret_tend           == host_data[s].nc_accret_tend);
      REQUIRE(back_to_cell_average_data[s].nc_selfcollect_tend      == host_data[s].nc_selfcollect_tend);
      REQUIRE(back_to_cell_average_data[s].nc2nr_autoconv_tend      == host_data[s].nc2nr_autoconv_tend);
      REQUIRE(back_to_cell_average_data[s].nr_selfcollect_tend      == host_data[s].nr_selfcollect_tend);
      REQUIRE(back_to_cell_average_data[s].nr_evap_tend             == host_data[s].nr_evap_tend);
      REQUIRE(back_to_cell_average_data[s].ncautr                   == host_data[s].ncautr);
      REQUIRE(back_to_cell_average_data[s].qi2qv_sublim_tend        == host_data[s].qi2qv_sublim_tend);
      REQUIRE(back_to_cell_average_data[s].nr_ice_shed_tend         == host_data[s].nr_ice_shed_tend);
      REQUIRE(back_to_cell_average_data[s].qc2qi_hetero_freeze_tend == host_data[s].qc2qi_hetero_freeze_tend);
      REQUIRE(back_to_cell_average_data[s].qr2qi_collect_tend       == host_data[s].qr2qi_collect_tend);
      REQUIRE(back_to_cell_average_data[s].qc2qr_ice_shed_tend      == host_data[s].qc2qr_ice_shed_tend);
      REQUIRE(back_to_cell_average_data[s].qi2qr_melt_tend          == host_data[s].qi2qr_melt_tend);
      REQUIRE(back_to_cell_average_data[s].qc2qi_collect_tend       == host_data[s].qc2qi_collect_tend);
      REQUIRE(back_to_cell_average_data[s].qr2qi_immers_freeze_tend == host_data[s].qr2qi_immers_freeze_tend);
      REQUIRE(back_to_cell_average_data[s].ni2nr_melt_tend          == host_data[s].ni2nr_melt_tend);
      REQUIRE(back_to_cell_average_data[s].nc_collect_tend          == host_data[s].nc_collect_tend);
      REQUIRE(back_to_cell_average_data[s].ncshdc                   == host_data[s].ncshdc);
      REQUIRE(back_to_cell_average_data[s].nc2ni_immers_freeze_tend == host_data[s].nc2ni_immers_freeze_tend);
      REQUIRE(back_to_cell_average_data[s].nr_collect_tend          == host_data[s].nr_collect_tend);
      REQUIRE(back_to_cell_average_data[s].ni_selfcollect_tend      == host_data[s].ni_selfcollect_tend);
      REQUIRE(back_to_cell_average_data[s].qv2qi_vapdep_tend        == host_data[s].qv2qi_vapdep_tend);
      REQUIRE(back_to_cell_average_data[s].nr2ni_immers_freeze_tend == host_data[s].nr2ni_immers_freeze_tend);
      REQUIRE(back_to_cell_average_data[s].ni_sublim_tend           == host_data[s].ni_sublim_tend);
      REQUIRE(back_to_cell_average_data[s].qv2qi_nucleat_tend       == host_data[s].qv2qi_nucleat_tend);
      REQUIRE(back_to_cell_average_data[s].ni_nucleat_tend          == host_data[s].ni_nucleat_tend);
      REQUIRE(back_to_cell_average_data[s].qc2qi_berg_tend          == host_data[s].qc2qi_berg_tend);
    }
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_back_to_cell_average", "[p3_functions]")
{
  using TRIF = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestBackToCellAverage;

  TRIF::run_phys();
  TRIF::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
