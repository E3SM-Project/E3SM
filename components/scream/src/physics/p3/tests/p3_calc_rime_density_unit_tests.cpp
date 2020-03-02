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
struct UnitWrap::UnitTest<D>::TestCalcRimeDensity {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);

  // This is the threshold for whether the qc cloud mixing ratio is large
  // enough to affect the rime density.
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar qc_incld_small = 0.9 * qsmall;
  constexpr Scalar qc_incld_not_small = 2.0 * qsmall;

  // We need to test the calculation under freezing and not-freezing
  // conditions.
  constexpr Scalar t_freezing = 0.9 * C::RainFrze,
                   t_not_freezing = 2.0 * C::RainFrze;

  // Ideally, we'd also test the calculation based on the mass-weighted mean
  // size Ri--whether it's above or below 8, specifically. Unfortunately,
  // computing Ri is complicated, and can't be done with black-box testing.
  // So maybe we run the gambit with various values given for other params
  // until we think of a better approach.
  constexpr Scalar rhofaci1 = 0.25, rhofaci2 = 0.5, rhofaci3 = 0.75, rhofaci4 = 1.0;
  constexpr Scalar f1pr021 = 0.125, f1pr022 = 0.25, f1pr023 = 0.375, f1pr024 = 0.5;
  constexpr Scalar acn1 = 0.1, acn2 = 0.4, acn3 = 0.8, acn4 = 0.9;
  constexpr Scalar lamc1 = 0.1, lamc2 = 0.2, lamc3 = 0.3, lamc4 = 0.4;
  constexpr Scalar mu_c1 = 0.2, mu_c2 = 0.4, mu_c3 = 0.6, mu_c4 = 0.8;
  constexpr Scalar qccol1 = 0.3, qccol2 = 0.4, qccol3 = 0.5, qccol4 = 0.6;

  CalcRimeDensityData calc_rime_density_data[max_pack_size] = {
    // t, lamc, mu_c, cdist1, qc_incld
    {t_not_freezing, rhofaci1, f1pr021, acn1, lamc1, mu_c1, qc_incld_small, qccol1},
    {t_not_freezing, rhofaci2, f1pr022, acn2, lamc2, mu_c2, qc_incld_small, qccol2},
    {t_not_freezing, rhofaci3, f1pr023, acn3, lamc3, mu_c3, qc_incld_small, qccol3},
    {t_not_freezing, rhofaci4, f1pr024, acn4, lamc4, mu_c4, qc_incld_small, qccol4},

    {t_not_freezing, rhofaci1, f1pr021, acn1, lamc1, mu_c1, qc_incld_not_small, qccol1},
    {t_not_freezing, rhofaci2, f1pr022, acn2, lamc2, mu_c2, qc_incld_not_small, qccol2},
    {t_not_freezing, rhofaci3, f1pr023, acn3, lamc3, mu_c3, qc_incld_not_small, qccol3},
    {t_not_freezing, rhofaci4, f1pr024, acn4, lamc4, mu_c4, qc_incld_not_small, qccol4},

    {t_freezing, rhofaci1, f1pr021, acn1, lamc1, mu_c1, qc_incld_small, qccol1},
    {t_freezing, rhofaci2, f1pr022, acn2, lamc2, mu_c2, qc_incld_small, qccol2},
    {t_freezing, rhofaci3, f1pr023, acn3, lamc3, mu_c3, qc_incld_small, qccol3},
    {t_freezing, rhofaci4, f1pr024, acn4, lamc4, mu_c4, qc_incld_small, qccol4},

    {t_freezing, rhofaci1, f1pr021, acn1, lamc1, mu_c1, qc_incld_not_small, qccol1},
    {t_freezing, rhofaci2, f1pr022, acn2, lamc2, mu_c2, qc_incld_not_small, qccol2},
    {t_freezing, rhofaci3, f1pr023, acn3, lamc3, mu_c3, qc_incld_not_small, qccol3},
    {t_freezing, rhofaci4, f1pr024, acn4, lamc4, mu_c4, qc_incld_not_small, qccol4}
  };

  // Sync to device
  view_1d<CalcRimeDensityData> device_data("calc_rime_density", Spack::n);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&calc_rime_density_data[0], &calc_rime_density_data[0] + Spack::n,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < Spack::n; ++i) {
    calc_rime_density(calc_rime_density_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack t, rhofaci, f1pr02, acn, lamc, mu_c, qc_incld, qccol;
    for (Int s = 0; s < Spack::n; ++s) {
      t[s]        = device_data(s).t;
      rhofaci[s]  = device_data(s).rhofaci;
      f1pr02[s]   = device_data(s).f1pr02;
      acn[s]      = device_data(s).acn;
      lamc[s]     = device_data(s).lamc;
      mu_c[s]     = device_data(s).mu_c;
      qc_incld[s] = device_data(s).qc_incld;
      qccol[s]    = device_data(s).qccol;
    }

    Spack vtrmi1{0.0};
    Spack rhorime_c{0.0};

    Functions::calc_rime_density(t, rhofaci, f1pr02, acn, lamc, mu_c,
                                 qc_incld, qccol, vtrmi1, rhorime_c);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      device_data(s).vtrmi1  = vtrmi1[s];
      device_data(s).rhorime_c  = rhorime_c[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < Spack::n; ++s) {
    REQUIRE(calc_rime_density_data[s].vtrmi1 == host_data[s].vtrmi1);
    REQUIRE(calc_rime_density_data[s].rhorime_c == host_data[s].rhorime_c);
  }

}

};

}
}
}

namespace {

TEST_CASE("p3_calc_rime_density", "[p3_functions]")
{
  using TRIF = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcRimeDensity;

  TRIF::run_phys();
  TRIF::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
