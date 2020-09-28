#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
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
struct UnitWrap::UnitTest<D>::TestPreventIceOverdepletion {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  constexpr Scalar p1 = 0.1, p2 = 0.2, p3 = 0.3, p4 = 0.4;
  constexpr Scalar t1 = 0.2, t2 = 0.4, t3 = 0.6, t4 = 0.8;
  constexpr Scalar qv1 = 0.125, qv2 = 0.25, qv3 = 0.375, qv4 = 0.5;
  constexpr Scalar latent_heat_sublim1 = 0.25, latent_heat_sublim2 = 0.5, latent_heat_sublim3 = 0.75, latent_heat_sublim4 = 1.0;
  constexpr Scalar inv_dt = 0.11;
  constexpr Scalar qv2qi_vapdep_tend = 1.0;
  constexpr Scalar qi2qv_sublim_tend = 1.0;

  PreventIceOverdepletionData prevent_ice_overdepletion_data[max_pack_size] = {
    // pres, T_atm, qv, latent_heat_sublim, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend
    {p1, t1, qv1, latent_heat_sublim1, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p1, t1, qv1, latent_heat_sublim1, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p1, t1, qv1, latent_heat_sublim1, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p1, t1, qv1, latent_heat_sublim1, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},

    {p2, t2, qv2, latent_heat_sublim2, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p2, t2, qv2, latent_heat_sublim2, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p2, t2, qv2, latent_heat_sublim2, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p2, t2, qv2, latent_heat_sublim2, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},

    {p3, t3, qv3, latent_heat_sublim3, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p3, t3, qv3, latent_heat_sublim3, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p3, t3, qv3, latent_heat_sublim3, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p3, t3, qv3, latent_heat_sublim3, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},

    {p4, t4, qv4, latent_heat_sublim4, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p4, t4, qv4, latent_heat_sublim4, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p4, t4, qv4, latent_heat_sublim4, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend},
    {p4, t4, qv4, latent_heat_sublim4, inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend}
  };

  // Sync to device
  view_1d<PreventIceOverdepletionData> device_data("prevent_ice_overdepletion", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&prevent_ice_overdepletion_data[0], &prevent_ice_overdepletion_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < max_pack_size; ++i) {
    prevent_ice_overdepletion(prevent_ice_overdepletion_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack pres, T_atm, qv, latent_heat_sublim, qv2qi_vapdep_tend, qi2qv_sublim_tend;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      pres[s]               = device_data(vs).pres;
      T_atm[s]                  = device_data(vs).T_atm;
      qv[s]                 = device_data(vs).qv;
      latent_heat_sublim[s] = device_data(vs).latent_heat_sublim;
      qv2qi_vapdep_tend[s]  = device_data(vs).qv2qi_vapdep_tend;
      qi2qv_sublim_tend[s]  = device_data(vs).qi2qv_sublim_tend;
    }

    Functions::prevent_ice_overdepletion(pres, T_atm, qv, latent_heat_sublim, device_data(0).inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).qv2qi_vapdep_tend  = qv2qi_vapdep_tend[s];
      device_data(vs).qi2qv_sublim_tend  = qi2qv_sublim_tend[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < max_pack_size; ++s) {
    REQUIRE(prevent_ice_overdepletion_data[s].qv2qi_vapdep_tend == host_data[s].qv2qi_vapdep_tend);
    REQUIRE(prevent_ice_overdepletion_data[s].qi2qv_sublim_tend == host_data[s].qi2qv_sublim_tend);
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_prevent_ice_overdepletion", "[p3_functions]")
{
  using TRIF = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPreventIceOverdepletion;

  TRIF::run_phys();
  TRIF::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
