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
struct UnitWrap::UnitTest<D>::TestIceSed {

static void run_phys_calc_bulk_rhime()
{
  // TODO
}

static void run_phys_ice_sed()
{
  // TODO
}

static void run_phys()
{
  run_phys_calc_bulk_rhime();
  run_phys_ice_sed();
}

static void run_bfb_calc_bulk_rhime()
{
  static constexpr Scalar qsmall = C::QSMALL;
  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);

  // Load some lookup inputs, need at least one per pack value
  CalcBulkRhoRimeData cbrr_fortran[max_pack_size] = {
    //     qi_tot,       qi_rim,       bi_rim
    {9.999978E-08, 9.999978E-03, 1.111108E-10},
    {0.000000E+00, 8.571428E-05, 1.000000E-02},
    {1.800685E-12, 1.818806E-13, 6.272458E-12},
    {5.164017E-10, 0.000000E+00, 0.000000E+00},

    {9.999978E-08, 0.000000E+00, 1.111108E-10},
    {5.100000E-03, 8.571428E-05, 1.000000E-02},
    {0.000000E+00, 1.818806E-13, 6.272458E-12},
    {5.164017E-10, 0.000000E+00, 0.000000E+00},

    {9.999978E-08, 9.999978E-08, 1.111108E-10},
    {5.100000E-03, 0.000000E+00, 1.000000E-02},
    {1.800685E-12, 1.818806E-13, 6.272458E-17},
    {0.000000E+00, 1.818806E-13, 0.000000E+00},

    {0.000000E+00, 9.999978E-08, 1.111108E-17},
    {5.100000E-03, 8.571428E-05, 1.000000E-02},
    {0.000000E+00, 1.818806E-13, 6.272458E-12},
    {5.164017E-10, 0.000000E+00, 0.000000E+00},
  };

  // Sync to device, needs to happen before fortran calls so that
  // inout data is in original state
  view_1d<CalcBulkRhoRimeData> cbrr_device("cbrr", Spack::n);
  const auto cbrr_host = Kokkos::create_mirror_view(cbrr_device);
  std::copy(&cbrr_fortran[0], &cbrr_fortran[0] + Spack::n, cbrr_host.data());
  Kokkos::deep_copy(cbrr_device, cbrr_host);

  // Get data from fortran
  for (Int i = 0; i < max_pack_size; ++i) {
    if (cbrr_fortran[i].qi_tot > qsmall) {
      calc_bulk_rho_rime(cbrr_fortran[i]);
    }
  }

  // Calc bulk rime from a kernel and copy results back to host
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    // Init pack inputs
    Spack qi_tot, qi_rim, bi_rim;
    for (Int s = 0; s < Spack::n; ++s) {
      qi_tot[s] = cbrr_device(s).qi_tot;
      qi_rim[s] = cbrr_device(s).qi_rim;
      bi_rim[s] = cbrr_device(s).bi_rim;
    }

    Smask gt_small(qi_tot > qsmall);
    Spack rho_rime = Functions::calc_bulk_rho_rime(gt_small, qi_tot, qi_rim, bi_rim);

    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      cbrr_device(s).qi_rim   = qi_rim[s];
      cbrr_device(s).bi_rim   = bi_rim[s];
      cbrr_device(s).rho_rime = rho_rime[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(cbrr_host, cbrr_device);

  // Validate results
  for (Int s = 0; s < Spack::n; ++s) {
    REQUIRE(cbrr_fortran[s].qi_rim   == cbrr_host(s).qi_rim);
    REQUIRE(cbrr_fortran[s].bi_rim   == cbrr_host(s).bi_rim);
    REQUIRE(cbrr_fortran[s].rho_rime == cbrr_host(s).rho_rime);
  }
}

static void run_bfb_ice_sed()
{
  // TODO
}

static void run_bfb()
{
  run_bfb_calc_bulk_rhime();
  run_bfb_ice_sed();
}

};

}
}
}

namespace {

TEST_CASE("p3_ice_sed", "[p3_functions]")
{
  using TCS = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceSed;

  TCS::run_phys();
  TCS::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
