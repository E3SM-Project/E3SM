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
struct UnitWrap::UnitTest<D>::TestCldliqImmersionFreezing {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
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
  constexpr Scalar qc_relvar_val = 1;

  CldliqImmersionFreezingData cldliq_imm_freezing_data[max_pack_size] = {
    // t, lamc, mu_c, cdist1, qc_incld, qc_relvar
    {t_not_freezing, lamc1, mu_c1, cdist11, qc_incld_small,qc_relvar_val},
    {t_not_freezing, lamc2, mu_c2, cdist12, qc_incld_small,qc_relvar_val},
    {t_not_freezing, lamc3, mu_c3, cdist13, qc_incld_small,qc_relvar_val},
    {t_not_freezing, lamc4, mu_c4, cdist14, qc_incld_small,qc_relvar_val},

    {t_not_freezing, lamc1, mu_c1, cdist11, qc_incld_not_small,qc_relvar_val},
    {t_not_freezing, lamc2, mu_c2, cdist12, qc_incld_not_small,qc_relvar_val},
    {t_not_freezing, lamc3, mu_c3, cdist13, qc_incld_not_small,qc_relvar_val},
    {t_not_freezing, lamc4, mu_c4, cdist14, qc_incld_not_small,qc_relvar_val},

    {t_freezing, lamc1, mu_c1, cdist11, qc_incld_small,qc_relvar_val},
    {t_freezing, lamc2, mu_c2, cdist12, qc_incld_small,qc_relvar_val},
    {t_freezing, lamc3, mu_c3, cdist13, qc_incld_small,qc_relvar_val},
    {t_freezing, lamc4, mu_c4, cdist14, qc_incld_small,qc_relvar_val},

    {t_freezing, lamc1, mu_c1, cdist11, qc_incld_not_small,qc_relvar_val},
    {t_freezing, lamc2, mu_c2, cdist12, qc_incld_not_small,qc_relvar_val},
    {t_freezing, lamc3, mu_c3, cdist13, qc_incld_not_small,qc_relvar_val},
    {t_freezing, lamc4, mu_c4, cdist14, qc_incld_not_small,qc_relvar_val}
  };

  // Sync to device
  view_1d<CldliqImmersionFreezingData> device_data("cldliq_imm_freezing", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&cldliq_imm_freezing_data[0], &cldliq_imm_freezing_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < max_pack_size; ++i) {
    cldliq_immersion_freezing(cldliq_imm_freezing_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack t, lamc, mu_c, cdist1, qc_incld,qc_relvar;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      t[s]        = device_data(vs).t;
      lamc[s]     = device_data(vs).lamc;
      mu_c[s]     = device_data(vs).mu_c;
      cdist1[s]   = device_data(vs).cdist1;
      qc_incld[s] = device_data(vs).qc_incld;
      qc_relvar[s]= device_data(vs).qc_relvar;
    }

    Spack qcheti{0.0};
    Spack ncheti{0.0};

    Functions::cldliq_immersion_freezing(t, lamc, mu_c, cdist1, qc_incld, qc_relvar,
                                         qcheti, ncheti);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).qcheti  = qcheti[s];
      device_data(vs).ncheti  = ncheti[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < max_pack_size; ++s) {
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
