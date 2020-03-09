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
struct UnitWrap::UnitTest<D>::TestPreventIceOverdepletion {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);

  constexpr Scalar p1 = 0.1, p2 = 0.2, p3 = 0.3, p4 = 0.4;
  constexpr Scalar t1 = 0.2, t2 = 0.4, t3 = 0.6, t4 = 0.8;
  constexpr Scalar qv1 = 0.125, qv2 = 0.25, qv3 = 0.375, qv4 = 0.5;
  constexpr Scalar xxls1 = 0.25, xxls2 = 0.5, xxls3 = 0.75, xxls4 = 1.0;
  constexpr Scalar odt1 = 0.11, odt2 = 0.22, odt3 = 0.33, odt4 = 0.44;

  PreventIceOverdepletionData prevent_ice_overdepletion_data[max_pack_size] = {
    // pres, t, qv, xxls, odt, qidep, qisub
    {p1, t1, qv1, xxls1, odt1},
    {p1, t1, qv1, xxls1, odt2},
    {p1, t1, qv1, xxls1, odt3},
    {p1, t1, qv1, xxls1, odt4},

    {p2, t2, qv2, xxls2, odt1},
    {p2, t2, qv2, xxls2, odt2},
    {p2, t2, qv2, xxls2, odt3},
    {p2, t2, qv2, xxls2, odt4},

    {p3, t3, qv3, xxls3, odt1},
    {p3, t3, qv3, xxls3, odt2},
    {p3, t3, qv3, xxls3, odt3},
    {p3, t3, qv3, xxls3, odt4},

    {p4, t4, qv4, xxls4, odt1},
    {p4, t4, qv4, xxls4, odt2},
    {p4, t4, qv4, xxls4, odt3},
    {p4, t4, qv4, xxls4, odt4}
  };

  // Sync to device
  view_1d<PreventIceOverdepletionData> device_data("prevent_ice_overdepletion", Spack::n);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&prevent_ice_overdepletion_data[0], &prevent_ice_overdepletion_data[0] + Spack::n,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Run the Fortran subroutine.
  for (Int i = 0; i < Spack::n; ++i) {
    prevent_ice_overdepletion(prevent_ice_overdepletion_data[i]);
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack pres, t, qv, xxls, odt;
    for (Int s = 0; s < Spack::n; ++s) {
      pres[s] = device_data(s).pres;
      t[s]    = device_data(s).t;
      qv[s]   = device_data(s).qv;
      xxls[s] = device_data(s).xxls;
      odt[s]  = device_data(s).odt;
    }

    Spack qidep{0.0};
    Spack qisub{0.0};

    Functions::prevent_ice_overdepletion(pres, t, qv, xxls, odt, qidep, qisub);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      device_data(s).qidep  = qidep[s];
      device_data(s).qisub  = qisub[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  for (Int s = 0; s < Spack::n; ++s) {
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
