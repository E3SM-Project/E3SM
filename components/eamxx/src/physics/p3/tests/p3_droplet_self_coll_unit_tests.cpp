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
struct UnitWrap::UnitTest<D>::TestDropletSelfCollection : public UnitWrap::UnitTest<D>::Base {

void run_phys()
{
  // TODO
}

void run_bfb()
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

  DropletSelfCollectionData droplet_self_coll_data[max_pack_size] = {
    // rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend, nc_selfcollect_tend
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

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < max_pack_size; ++i) {
      droplet_self_coll_data[i].read(Base::m_fid);
    }
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      rho[s]                   = device_data(vs).rho;
      inv_rho[s]               = device_data(vs).inv_rho;
      qc_incld[s]              = device_data(vs).qc_incld;
      mu_c[s]                  = device_data(vs).mu_c;
      nu[s]                    = device_data(vs).nu;
      nc2nr_autoconv_tend[s]   = device_data(vs).nc2nr_autoconv_tend;
    }

    Spack nc_selfcollect_tend{0.0};

    Functions::droplet_self_collection(rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend,
                                       nc_selfcollect_tend);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).nc_selfcollect_tend  = nc_selfcollect_tend[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(droplet_self_coll_data[s].nc_selfcollect_tend == host_data[s].nc_selfcollect_tend);
    }
  }
  else if (this->m_baseline_action == GENERATE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      host_data(s).write(Base::m_fid);
    }
  }

}

};

}
}
}

namespace {

TEST_CASE("p3_droplet_self_collection", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDropletSelfCollection;

  T t;
  t.run_phys();
  t.run_bfb();
}

} // namespace
