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
struct UnitWrap::UnitTest<D>::TestCalcRimeDensity : public UnitWrap::UnitTest<D>::Base {

void run_phys()
{
  // TODO
}

void run_bfb()
{
  // This is the threshold for whether the qc cloud mixing ratio is large
  // enough to affect the rime density.
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar qc_incld_small = 0.9 * qsmall;
  constexpr Scalar qc_incld_not_small = 2.0 * qsmall;

  // Same goes for qc2qi_collect_tend.
  constexpr Scalar qc2qi_collect_tend_small = 0.9 * qsmall;
  constexpr Scalar qc2qi_collect_tend_not_small = 2.0 * qsmall;

  // We need to test the calculation under freezing and not-freezing
  // conditions.
  constexpr Scalar t_freezing = 0.9 * C::T_rainfrz,
                   t_not_freezing = 2.0 * C::T_rainfrz;

  // Ideally, we'd also test the calculation based on the mass-weighted mean
  // size Ri--whether it's above or below 8, specifically. Unfortunately,
  // computing Ri is complicated, and can't be done with black-box testing.
  // So maybe we run the gambit with various values given for other params
  // until we think of a better approach.
  constexpr Scalar rhofaci1 = 0.25, rhofaci2 = 0.5;
  constexpr Scalar table_val_qi_fallspd1 = 0.125, table_val_qi_fallspd2 = 0.875;
  constexpr Scalar acn1 = 0.1, acn2 = 0.6;
  constexpr Scalar lamc1 = 0.4, lamc2 = 0.8;
  constexpr Scalar mu_c1 = 0.2, mu_c2 = 0.7;

  CalcRimeDensityData calc_rime_density_data[max_pack_size] = {
    // T_atm, rhofaci, table_val_qi_fallspd1, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend
    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_small},

    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_not_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_not_small},

    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_small},

    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_not_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_not_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_not_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_not_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_not_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_not_small},
  };

  // Sync to device
  view_1d<CalcRimeDensityData> device_data("calc_rime_density", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&calc_rime_density_data[0], &calc_rime_density_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < max_pack_size; ++i) {
      calc_rime_density_data[i].read(Base::m_fid);
    }
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      T_atm[s]                = device_data(vs).T_atm;
      rhofaci[s]              = device_data(vs).rhofaci;
      table_val_qi_fallspd[s] = device_data(vs).table_val_qi_fallspd;
      acn[s]                  = device_data(vs).acn;
      lamc[s]                 = device_data(vs).lamc;
      mu_c[s]                 = device_data(vs).mu_c;
      qc_incld[s]             = device_data(vs).qc_incld;
      qc2qi_collect_tend[s]   = device_data(vs).qc2qi_collect_tend;
    }

    Spack vtrmi1{0.0};
    Spack rho_qm_cloud{0.0};

    Functions::calc_rime_density(T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c,
                                 qc_incld, qc2qi_collect_tend, vtrmi1, rho_qm_cloud);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).vtrmi1  = vtrmi1[s];
      device_data(vs).rho_qm_cloud  = rho_qm_cloud[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(calc_rime_density_data[s].vtrmi1 == host_data[s].vtrmi1);
      REQUIRE(calc_rime_density_data[s].rho_qm_cloud == host_data[s].rho_qm_cloud);
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

TEST_CASE("p3_calc_rime_density", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcRimeDensity;

  T t;
  t.run_phys();
  t.run_bfb();
}

} // namespace
