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
struct UnitWrap::UnitTest<D>::TestCldliqImmersionFreezing {

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
  constexpr Scalar qc_incld_small = 0.9 * qsmall;
  constexpr Scalar qc_incld_not_small = 2.0 * qsmall;
  constexpr Scalar lamc1 = 0.1, lamc2 = 0.2, lamc3 = 0.3, lamc4 = 0.4;
  constexpr Scalar mu_c1 = 0.2, mu_c2 = 0.4, mu_c3 = 0.6, mu_c4 = 0.8;
  constexpr Scalar cdist11 = 0.25, cdist12 = 0.5, cdist13 = 0.75, cdist14 = 1.0;

  CldliqImmersionFreezingData cldliq_imm_freezing_data[max_pack_size] = {
    // t, lamc, mu_c, cdist1, qc_incld
    {t_not_freezing, lamc1, mu_c1, cdist11, qc_incld_small},
    {t_not_freezing, lamc2, mu_c2, cdist12, qc_incld_small},
    {t_not_freezing, lamc3, mu_c3, cdist13, qc_incld_small},
    {t_not_freezing, lamc4, mu_c4, cdist14, qc_incld_small},

    {t_not_freezing, lamc1, mu_c1, cdist11, qc_incld_not_small},
    {t_not_freezing, lamc2, mu_c2, cdist12, qc_incld_not_small},
    {t_not_freezing, lamc3, mu_c3, cdist13, qc_incld_not_small},
    {t_not_freezing, lamc4, mu_c4, cdist14, qc_incld_not_small},

    {t_freezing, lamc1, mu_c1, cdist11, qc_incld_small},
    {t_freezing, lamc2, mu_c2, cdist12, qc_incld_small},
    {t_freezing, lamc3, mu_c3, cdist13, qc_incld_small},
    {t_freezing, lamc4, mu_c4, cdist14, qc_incld_small},

    {t_freezing, lamc1, mu_c1, cdist11, qc_incld_not_small},
    {t_freezing, lamc2, mu_c2, cdist12, qc_incld_not_small},
    {t_freezing, lamc3, mu_c3, cdist13, qc_incld_not_small},
    {t_freezing, lamc4, mu_c4, cdist14, qc_incld_not_small}
  };

  // Sync to device
  view_1d<CldliqImmersionFreezingData> device_data("cldliq_imm_freezing", Spack::n);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&cldliq_imm_freezing_data[0], &cldliq_imm_freezing_data[0] + Spack::n,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < Spack::n; ++i) {
    cldliq_immersion_freezing(cldliq_imm_freezing_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack t, lamc, mu_c, cdist1, qc_incld;
    for (Int s = 0; s < Spack::n; ++s) {
      t[s]        = device_data(s).t;
      lamc[s]     = device_data(s).lamc;
      mu_c[s]     = device_data(s).mu_c;
      cdist1[s]   = device_data(s).cdist1;
      qc_incld[s] = device_data(s).qc_incld;
    }

    Spack qcheti{0.0};
    Spack ncheti{0.0};

    Functions::cldliq_immersion_freezing(t, lamc, mu_c, cdist1, qc_incld,
                                         qcheti, ncheti);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      device_data(s).qcheti  = qcheti[s];
      device_data(s).ncheti  = ncheti[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < Spack::n; ++s) {
    REQUIRE(cldliq_imm_freezing_data[s].qcheti == host_data[s].qcheti);
    REQUIRE(cldliq_imm_freezing_data[s].ncheti == host_data[s].ncheti);
  }

}

};

}
}
}

namespace {

TEST_CASE("p3_cldliq_immersion_freezing", "[p3_functions]")
{
  using TRIF = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCldliqImmersionFreezing;

  TRIF::run_phys();
  TRIF::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
