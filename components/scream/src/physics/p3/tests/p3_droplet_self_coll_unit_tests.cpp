#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
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
struct UnitWrap::UnitTest<D>::TestDropletSelfCollection {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  // This is the threshold for whether the qc and qr cloud mixing ratios are
  // large enough to affect the warm-phase process rates qcacc and ncacc.
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

  DropletSelfCollectionData droplet_self_coll_data[max_pack_size] = {
    // rho, inv_rho, qc_incld, mu_c, nu, ncautc, ncslf
    {rho1, inv_rho1, qc_incld_small, nc_incld1, qr_incld_small},
    {rho2, inv_rho2, qc_incld_small, nc_incld2, qr_incld_small},
    {rho3, inv_rho3, qc_incld_small, nc_incld3, qr_incld_small},
    {rho4, inv_rho4, qc_incld_small, nc_incld4, qr_incld_small},

    {rho1, inv_rho1, qc_incld_small, nc_incld1, qr_incld_not_small},
    {rho2, inv_rho2, qc_incld_small, nc_incld2, qr_incld_not_small},
    {rho3, inv_rho3, qc_incld_small, nc_incld3, qr_incld_not_small},
    {rho4, inv_rho4, qc_incld_small, nc_incld4, qr_incld_not_small},

    {rho1, inv_rho1, qc_incld_not_small, nc_incld1, qr_incld_small},
    {rho2, inv_rho2, qc_incld_not_small, nc_incld2, qr_incld_small},
    {rho3, inv_rho3, qc_incld_not_small, nc_incld3, qr_incld_small},
    {rho4, inv_rho4, qc_incld_not_small, nc_incld4, qr_incld_small},

    {rho1, inv_rho1, qc_incld_not_small, nc_incld1, qr_incld_not_small},
    {rho2, inv_rho2, qc_incld_not_small, nc_incld2, qr_incld_not_small},
    {rho3, inv_rho3, qc_incld_not_small, nc_incld3, qr_incld_not_small},
    {rho4, inv_rho4, qc_incld_not_small, nc_incld4, qr_incld_not_small}
  };

  // Sync to device
  view_1d<DropletSelfCollectionData> device_data("droplet_self_coll", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&droplet_self_coll_data[0], &droplet_self_coll_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < max_pack_size; ++i) {
    droplet_self_collection(droplet_self_coll_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack rho, inv_rho, qc_incld, mu_c, nu, ncautc;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      rho[s]      = device_data(vs).rho;
      inv_rho[s]  = device_data(vs).inv_rho;
      qc_incld[s] = device_data(vs).qc_incld;
      mu_c[s]     = device_data(vs).mu_c;
      nu[s]       = device_data(vs).nu;
      ncautc[s]   = device_data(vs).ncautc;
    }

    Spack ncslf{0.0};

    Functions::droplet_self_collection(rho, inv_rho, qc_incld, mu_c, nu, ncautc,
                                       ncslf);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).ncslf  = ncslf[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < max_pack_size; ++s) {
    REQUIRE(droplet_self_coll_data[s].ncslf == host_data[s].ncslf);
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_droplet_self_collection", "[p3_functions]")
{
  using TCRA = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDropletSelfCollection;

  TCRA::run_phys();
  TCRA::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
