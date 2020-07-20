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
  constexpr Scalar xxls1 = 0.25, xxls2 = 0.5, xxls3 = 0.75, xxls4 = 1.0;
  constexpr Scalar odt = 0.11;
  constexpr Scalar qidep = 1.0;
  constexpr Scalar qisub = 1.0;

  PreventIceOverdepletionData prevent_ice_overdepletion_data[max_pack_size] = {
    // pres, t, qv, xxls, odt, qidep, qisub
    {p1, t1, qv1, xxls1, odt, qidep, qisub},
    {p1, t1, qv1, xxls1, odt, qidep, qisub},
    {p1, t1, qv1, xxls1, odt, qidep, qisub},
    {p1, t1, qv1, xxls1, odt, qidep, qisub},

    {p2, t2, qv2, xxls2, odt, qidep, qisub},
    {p2, t2, qv2, xxls2, odt, qidep, qisub},
    {p2, t2, qv2, xxls2, odt, qidep, qisub},
    {p2, t2, qv2, xxls2, odt, qidep, qisub},

    {p3, t3, qv3, xxls3, odt, qidep, qisub},
    {p3, t3, qv3, xxls3, odt, qidep, qisub},
    {p3, t3, qv3, xxls3, odt, qidep, qisub},
    {p3, t3, qv3, xxls3, odt, qidep, qisub},

    {p4, t4, qv4, xxls4, odt, qidep, qisub},
    {p4, t4, qv4, xxls4, odt, qidep, qisub},
    {p4, t4, qv4, xxls4, odt, qidep, qisub},
    {p4, t4, qv4, xxls4, odt, qidep, qisub}
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
    Spack pres, t, qv, xxls, qidep, qisub;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      pres[s]  = device_data(vs).pres;
      t[s]     = device_data(vs).t;
      qv[s]    = device_data(vs).qv;
      xxls[s]  = device_data(vs).xxls;
      qidep[s] = device_data(vs).qidep;
      qisub[s] = device_data(vs).qisub;
    }

    Functions::prevent_ice_overdepletion(pres, t, qv, xxls, device_data(0).odt, qidep, qisub);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      device_data(vs).qidep  = qidep[s];
      device_data(vs).qisub  = qisub[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < max_pack_size; ++s) {
    REQUIRE(prevent_ice_overdepletion_data[s].qidep == host_data[s].qidep);
    REQUIRE(prevent_ice_overdepletion_data[s].qisub == host_data[s].qisub);
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
