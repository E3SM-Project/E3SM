#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestCloudRainAccretion {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  // This is the threshold for whether the qc and qr cloud mixing ratios are
  // large enough to affect the warm-phase process rates qc2qr_accret_tend and nc_accret_tend.
  constexpr Scalar qsmall = C::QSMALL;

  constexpr Scalar rho1 = 4.056E-03, rho2 = 6.852E-02,
                   rho3 = 8.852E-02, rho4 = 1.902E-01;
  constexpr Scalar inv_rho1 = 1.0/rho1, inv_rho2 = 1.0/rho2,
                   inv_rho3 = 1.0/rho3, inv_rho4 = 1.0/rho4;
  constexpr Scalar qc_incld_small = 0.9 * qsmall;
  constexpr Scalar qr_incld_small = 0.9 * qsmall;
  constexpr Scalar qc_incld_not_small = 2.0 * qsmall;
  constexpr Scalar qr_incld_not_small = 2.0 * qsmall;
  constexpr Scalar nc_incld1 = 9.952E+05, nc_incld2 = 9.952E+06,
                   nc_incld3 = 1.734E+07, nc_incld4 = 9.952E+08;
  constexpr Scalar inv_qc_relvar_val = 1;

  CloudRainAccretionData cloud_rain_acc_data[max_pack_size] = {
    // rho, inv_rho, qc_incld, nc_incld, qr_incld, qc2qr_accret_tend, nc_accret_tend, inv_qc_relvar
    {rho1, inv_rho1, qc_incld_small, nc_incld1, qr_incld_small,inv_qc_relvar_val},
    {rho2, inv_rho2, qc_incld_small, nc_incld2, qr_incld_small,inv_qc_relvar_val},
    {rho3, inv_rho3, qc_incld_small, nc_incld3, qr_incld_small,inv_qc_relvar_val},
    {rho4, inv_rho4, qc_incld_small, nc_incld4, qr_incld_small,inv_qc_relvar_val},

    {rho1, inv_rho1, qc_incld_small, nc_incld1, qr_incld_not_small,inv_qc_relvar_val},
    {rho2, inv_rho2, qc_incld_small, nc_incld2, qr_incld_not_small,inv_qc_relvar_val},
    {rho3, inv_rho3, qc_incld_small, nc_incld3, qr_incld_not_small,inv_qc_relvar_val},
    {rho4, inv_rho4, qc_incld_small, nc_incld4, qr_incld_not_small,inv_qc_relvar_val},

    {rho1, inv_rho1, qc_incld_not_small, nc_incld1, qr_incld_small,inv_qc_relvar_val},
    {rho2, inv_rho2, qc_incld_not_small, nc_incld2, qr_incld_small,inv_qc_relvar_val},
    {rho3, inv_rho3, qc_incld_not_small, nc_incld3, qr_incld_small,inv_qc_relvar_val},
    {rho4, inv_rho4, qc_incld_not_small, nc_incld4, qr_incld_small,inv_qc_relvar_val},

    {rho1, inv_rho1, qc_incld_not_small, nc_incld1, qr_incld_not_small,inv_qc_relvar_val},
    {rho2, inv_rho2, qc_incld_not_small, nc_incld2, qr_incld_not_small,inv_qc_relvar_val},
    {rho3, inv_rho3, qc_incld_not_small, nc_incld3, qr_incld_not_small,inv_qc_relvar_val},
    {rho4, inv_rho4, qc_incld_not_small, nc_incld4, qr_incld_not_small,inv_qc_relvar_val}
  };

  // Sync to device
  view_1d<CloudRainAccretionData> device_data("cloud_rain_acc", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&cloud_rain_acc_data[0], &cloud_rain_acc_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < max_pack_size; ++i) {
    cloud_rain_accretion(cloud_rain_acc_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack rho, inv_rho, qc_incld, nc_incld, qr_incld, inv_qc_relvar;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      rho[s]            = device_data(vs).rho;
      inv_rho[s]        = device_data(vs).inv_rho;
      qc_incld[s]       = device_data(vs).qc_incld;
      nc_incld[s]       = device_data(vs).nc_incld;
      qr_incld[s]       = device_data(vs).qr_incld;
      inv_qc_relvar[s]  = device_data(vs).inv_qc_relvar;
    }

    Spack qc2qr_accret_tend{0.0};
    Spack nc_accret_tend{0.0};

    Functions::cloud_rain_accretion(rho, inv_rho, qc_incld, nc_incld, qr_incld,
                                    inv_qc_relvar, qc2qr_accret_tend, nc_accret_tend);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).qc2qr_accret_tend  = qc2qr_accret_tend[s];
      device_data(vs).nc_accret_tend     = nc_accret_tend[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  if (SCREAM_BFB_TESTING) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(cloud_rain_acc_data[s].qc2qr_accret_tend == host_data[s].qc2qr_accret_tend);
      REQUIRE(cloud_rain_acc_data[s].nc_accret_tend == host_data[s].nc_accret_tend);
    }
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_cloud_rain_accretion", "[p3_functions]")
{
  using TCRA = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCloudRainAccretion;

  TCRA::run_phys();
  TCRA::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
