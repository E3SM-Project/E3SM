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
struct UnitWrap::UnitTest<D>::TestRainSed {

static void run_phys_rain_vel()
{
  // TODO
}

static void run_phys_rain_sed()
{
  // TODO
}

static void run_phys()
{
  run_phys_rain_vel();
  run_phys_rain_sed();
}

static void run_bfb_rain_vel()
{
  // Read in tables
  view_2d_table vn_table; view_2d_table vm_table; view_1d_table mu_r_table; view_dnu_table dnu;
  Functions::init_kokkos_tables(vn_table, vm_table, mu_r_table, dnu);

  constexpr Scalar qsmall = C::QSMALL;
  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);

  // Load some lookup inputs, need at least one per pack value
  ComputeRainFallVelocityData crfv_fortran[max_pack_size] = {
    // qr_incld,       rcldm,    rhofacr,         nr, nr_incld
    {1.1030E-04, 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {5.6298E-05, 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {1.0000E-02, 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},

    {0.0       , 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {0.0       , 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {1.0000E-02, 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},

    {1.1030E-04, 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {0.0       , 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {0.0       , 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},

    {0.0       , 1.0000E+00, 1.3221E+00, 6.5282E+07, 6.2964E+05},
    {2.1437E-13, 2.0000E+00, 1.0918E+00, 1.7935E+02, 6.5337E+07},
    {5.6298E-05, 3.0000E+00, 1.1129E+00, 6.2909E+05, 1.6576E+02},
    {1.0000E-02, 4.0000E+00, 1.0774E+00, 1.6986E+02, 1.9436E+02},
  };

  // Sync to device, needs to happen before fortran calls so that
  // inout data is in original state
  view_1d<ComputeRainFallVelocityData> crfv_device("crfv", Spack::n);
  const auto crfv_host = Kokkos::create_mirror_view(crfv_device);
  std::copy(&crfv_fortran[0], &crfv_fortran[0] + Spack::n, crfv_host.data());
  Kokkos::deep_copy(crfv_device, crfv_host);

  // Get data from fortran
  for (Int i = 0; i < max_pack_size; ++i) {
    if (crfv_fortran[i].qr_incld > qsmall) {
      compute_rain_fall_velocity(crfv_fortran[i]);
    }
  }

  // Calc bulk rime from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack qr_incld, rcldm, rhofacr, nr, nr_incld;
    for (Int s = 0; s < Spack::n; ++s) {
      qr_incld[s] = crfv_device(s).qr_incld;
      rcldm[s]    = crfv_device(s).rcldm;
      rhofacr[s]  = crfv_device(s).rhofacr;
      nr[s]       = crfv_device(s).nr;
      nr_incld[s] = crfv_device(s).nr_incld;
    }

    Smask gt_small(qr_incld > qsmall);
    Spack mu_r, lamr, V_qr, V_nr;
    Functions::compute_rain_fall_velocity(
      gt_small, vn_table, vm_table, qr_incld, rcldm, rhofacr, nr, nr_incld, mu_r, lamr, V_qr, V_nr);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      crfv_device(s).nr       = nr[s];
      crfv_device(s).nr_incld = nr_incld[s];
      crfv_device(s).mu_r     = mu_r[s];
      crfv_device(s).lamr     = lamr[s];
      crfv_device(s).V_qr     = V_qr[s];
      crfv_device(s).V_nr     = V_nr[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(crfv_host, crfv_device);

  // Validate results
  for (Int s = 0; s < Spack::n; ++s) {
    REQUIRE(crfv_fortran[s].nr       == crfv_host(s).nr);
    REQUIRE(crfv_fortran[s].nr_incld == crfv_host(s).nr_incld);
    REQUIRE(crfv_fortran[s].mu_r     == crfv_host(s).mu_r);
    REQUIRE(crfv_fortran[s].lamr     == crfv_host(s).lamr);
    REQUIRE(crfv_fortran[s].V_qr     == crfv_host(s).V_qr);
    REQUIRE(crfv_fortran[s].V_nr     == crfv_host(s).V_nr);
  }
}

static void run_bfb_rain_sed()
{
  // TODO
}

static void run_bfb()
{
  run_bfb_rain_vel();
  run_bfb_rain_sed();
}

};

}
}
}

namespace {

TEST_CASE("p3_rain_sed", "[p3_functions]")
{
  using TRS = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRainSed;

  TRS::run_phys();
  TRS::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
