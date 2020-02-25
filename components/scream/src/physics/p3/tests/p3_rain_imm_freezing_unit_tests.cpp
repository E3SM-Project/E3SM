#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestRainImmersionFreezing {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);

  // This is the threshold for whether the qc and qr cloud mixing ratios are
  // large enough to affect the warm-phase process rates qcacc and ncacc.
  constexpr Scalar qsmall = C::QSMALL;

  constexpr Scalar t_freezing = 0.9 * C::RainFrze,
                   t_not_freezing = 2.0 * C::RainFrze;
  constexpr Scalar qr_incld_small = 0.9 * qsmall;
  constexpr Scalar qr_incld_not_small = 2.0 * qsmall;
  constexpr Scalar lamr1 = 0.1, lamr2 = 0.2, lamr3 = 0.3, lamr4 = 0.4;
  constexpr Scalar mu_r1 = 0.2, mu_r2 = 0.4, mu_r3 = 0.6, mu_r4 = 0.8;
  constexpr Scalar cdistr1 = 0.25, cdistr2 = 0.5, cdistr3 = 0.75, cdistr4 = 1.0;

  RainImmersionFreezingData rain_imm_freezing_data[max_pack_size] = {
    // t, lamr, mu_r, cdistr, qr_incld, qrheti, nrheti
    {t_not_freezing, lamr1, mu_r1, cdistr1, qr_incld_small, 0.0, 0.0},
    {t_not_freezing, lamr2, mu_r2, cdistr2, qr_incld_small, 0.0, 0.0},
    {t_not_freezing, lamr3, mu_r3, cdistr3, qr_incld_small, 0.0, 0.0},
    {t_not_freezing, lamr4, mu_r4, cdistr4, qr_incld_small, 0.0, 0.0},

    {t_not_freezing, lamr1, mu_r1, cdistr1, qr_incld_not_small, 0.0, 0.0},
    {t_not_freezing, lamr2, mu_r2, cdistr2, qr_incld_not_small, 0.0, 0.0},
    {t_not_freezing, lamr3, mu_r3, cdistr3, qr_incld_not_small, 0.0, 0.0},
    {t_not_freezing, lamr4, mu_r4, cdistr4, qr_incld_not_small, 0.0, 0.0},

    {t_freezing, lamr1, mu_r1, cdistr1, qr_incld_small, 0.0, 0.0},
    {t_freezing, lamr2, mu_r2, cdistr2, qr_incld_small, 0.0, 0.0},
    {t_freezing, lamr3, mu_r3, cdistr3, qr_incld_small, 0.0, 0.0},
    {t_freezing, lamr4, mu_r4, cdistr4, qr_incld_small, 0.0, 0.0},

    {t_freezing, lamr1, mu_r1, cdistr1, qr_incld_not_small, 0.0, 0.0},
    {t_freezing, lamr2, mu_r2, cdistr2, qr_incld_not_small, 0.0, 0.0},
    {t_freezing, lamr3, mu_r3, cdistr3, qr_incld_not_small, 0.0, 0.0},
    {t_freezing, lamr4, mu_r4, cdistr4, qr_incld_not_small, 0.0, 0.0}
  };

  // Sync to device
  view_1d<RainImmersionFreezingData> device_data("rain_imm_freezing", Spack::n);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&rain_imm_freezing_data[0], &rain_imm_freezing_data[0] + Spack::n,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < Spack::n; ++i) {
    rain_immersion_freezing(rain_imm_freezing_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack t, lamr, mu_r, cdistr, qr_incld;
    for (Int s = 0; s < Spack::n; ++s) {
      t[s]        = device_data(s).t;
      lamr[s]     = device_data(s).lamr;
      mu_r[s]     = device_data(s).mu_r;
      cdistr[s]   = device_data(s).cdistr;
      qr_incld[s] = device_data(s).qr_incld;
    }

    Spack qrheti{0.0};
    Spack nrheti{0.0};

    Functions::rain_immersion_freezing(t, lamr, mu_r, cdistr, qr_incld,
                                       qrheti, nrheti);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      device_data(s).qrheti  = qrheti[s];
      device_data(s).nrheti  = nrheti[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < Spack::n; ++s) {
    REQUIRE(rain_imm_freezing_data[s].qrheti == host_data[s].qrheti);
    REQUIRE(rain_imm_freezing_data[s].nrheti == host_data[s].nrheti);
  }

}

};

}
}
}

namespace {

TEST_CASE("p3_cloud_rain_accretion", "[p3_functions]")
{
  using TRIF = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRainImmersionFreezing;

  TRIF::run_phys();
  TRIF::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
