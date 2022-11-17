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
struct UnitWrap::UnitTest<D>::TestRainImmersionFreezing {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  // This is the threshold for whether the qc and qr cloud mixing ratios are
  // large enough to affect the warm-phase process rates qc2qr_accret_tend and nc_accret_tend.
  constexpr Scalar qsmall = C::QSMALL;

  constexpr Scalar t_freezing = 0.9 * C::T_rainfrz,
                   t_not_freezing = 2.0 * C::T_rainfrz;
  constexpr Scalar qr_incld_small = 0.9 * qsmall;
  constexpr Scalar qr_incld_not_small = 2.0 * qsmall;
  constexpr Scalar lamr1 = 0.1, lamr2 = 0.2, lamr3 = 0.3, lamr4 = 0.4;
  constexpr Scalar mu_r1 = 0.2, mu_r2 = 0.4, mu_r3 = 0.6, mu_r4 = 0.8;
  constexpr Scalar cdistr1 = 0.25, cdistr2 = 0.5, cdistr3 = 0.75, cdistr4 = 1.0;

  RainImmersionFreezingData rain_imm_freezing_data[max_pack_size] = {
    // T_atm, lamr, mu_r, cdistr, qr_incld, qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend
    {t_not_freezing, lamr1, mu_r1, cdistr1, qr_incld_small},
    {t_not_freezing, lamr2, mu_r2, cdistr2, qr_incld_small},
    {t_not_freezing, lamr3, mu_r3, cdistr3, qr_incld_small},
    {t_not_freezing, lamr4, mu_r4, cdistr4, qr_incld_small},

    {t_not_freezing, lamr1, mu_r1, cdistr1, qr_incld_not_small},
    {t_not_freezing, lamr2, mu_r2, cdistr2, qr_incld_not_small},
    {t_not_freezing, lamr3, mu_r3, cdistr3, qr_incld_not_small},
    {t_not_freezing, lamr4, mu_r4, cdistr4, qr_incld_not_small},

    {t_freezing, lamr1, mu_r1, cdistr1, qr_incld_small},
    {t_freezing, lamr2, mu_r2, cdistr2, qr_incld_small},
    {t_freezing, lamr3, mu_r3, cdistr3, qr_incld_small},
    {t_freezing, lamr4, mu_r4, cdistr4, qr_incld_small},

    {t_freezing, lamr1, mu_r1, cdistr1, qr_incld_not_small},
    {t_freezing, lamr2, mu_r2, cdistr2, qr_incld_not_small},
    {t_freezing, lamr3, mu_r3, cdistr3, qr_incld_not_small},
    {t_freezing, lamr4, mu_r4, cdistr4, qr_incld_not_small}
  };

  // Sync to device
  view_1d<RainImmersionFreezingData> device_data("rain_imm_freezing", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&rain_imm_freezing_data[0], &rain_imm_freezing_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < max_pack_size; ++i) {
    rain_immersion_freezing(rain_imm_freezing_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack T_atm, lamr, mu_r, cdistr, qr_incld;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      T_atm[s]    = device_data(vs).T_atm;
      lamr[s]     = device_data(vs).lamr;
      mu_r[s]     = device_data(vs).mu_r;
      cdistr[s]   = device_data(vs).cdistr;
      qr_incld[s] = device_data(vs).qr_incld;
    }

    Spack qr2qi_immers_freeze_tend{0.0};
    Spack nr2ni_immers_freeze_tend{0.0};

    Functions::rain_immersion_freezing(T_atm, lamr, mu_r, cdistr, qr_incld,
                                       qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).qr2qi_immers_freeze_tend  = qr2qi_immers_freeze_tend[s];
      device_data(vs).nr2ni_immers_freeze_tend  = nr2ni_immers_freeze_tend[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  if (SCREAM_BFB_TESTING) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(rain_imm_freezing_data[s].qr2qi_immers_freeze_tend == host_data[s].qr2qi_immers_freeze_tend);
      REQUIRE(rain_imm_freezing_data[s].nr2ni_immers_freeze_tend == host_data[s].nr2ni_immers_freeze_tend);
    }
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_rain_immersion_freezing", "[p3_functions]")
{
  using TRIF = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRainImmersionFreezing;

  TRIF::run_phys();
  TRIF::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
