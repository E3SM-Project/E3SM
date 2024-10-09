#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

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
  view_2d_table vn_table_vals; view_2d_table vm_table_vals; view_2d_table revap_table_vals;
  view_1d_table mu_r_table_vals; view_dnu_table dnu;
  Functions::init_kokkos_tables(vn_table_vals, vm_table_vals, revap_table_vals, mu_r_table_vals, dnu);

  // Load some lookup inputs, need at least one per pack value
  ComputeRainFallVelocityData crfv_fortran[max_pack_size] = {
    // qr_incld,   rhofacr,   nr_incld
    {1.1030E-04, 1.3221E+00, 6.2964E+05},
    {2.1437E-13, 1.0918E+00, 6.5337E+07},
    {5.6298E-05, 1.1129E+00, 1.6576E+02},
    {1.0000E-02, 1.0774E+00, 1.9436E+02},

    {1.1030E-04, 1.3221E+00, 6.2964E+05},
    {2.1437E-13, 1.0918E+00, 6.5337E+07},
    {5.6298E-05, 1.1129E+00, 1.6576E+02},
    {1.0000E-02, 1.0774E+00, 1.9436E+02},

    {1.1030E-04, 1.3221E+00, 6.2964E+05},
    {2.1437E-13, 1.0918E+00, 6.5337E+07},
    {5.6298E-05, 1.1129E+00, 1.6576E+02},
    {1.0000E-02, 1.0774E+00, 1.9436E+02},

    {1.1030E-04, 1.3221E+00, 6.2964E+05},
    {2.1437E-13, 1.0918E+00, 6.5337E+07},
    {5.6298E-05, 1.1129E+00, 1.6576E+02},
    {1.0000E-02, 1.0774E+00, 1.9436E+02},

  };

  // Sync to device, needs to happen before fortran calls so that
  // inout data is in original state
  view_1d<ComputeRainFallVelocityData> crfv_device("crfv", max_pack_size);
  const auto crfv_host = Kokkos::create_mirror_view(crfv_device);
  std::copy(&crfv_fortran[0], &crfv_fortran[0] + max_pack_size, crfv_host.data());
  Kokkos::deep_copy(crfv_device, crfv_host);

  // Get data from fortran
  for (Int i = 0; i < max_pack_size; ++i) {
    compute_rain_fall_velocity(crfv_fortran[i]);
  }

  // Calc bulk rime from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack qr_incld, cld_frac_r, rhofacr, nr_incld;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      qr_incld[s] = crfv_device(vs).qr_incld;
      rhofacr[s]  = crfv_device(vs).rhofacr;
      nr_incld[s] = crfv_device(vs).nr_incld;
    }

    Spack mu_r(0), lamr(0), V_qr(0), V_nr(0);
    Functions::compute_rain_fall_velocity(
        vn_table_vals, vm_table_vals, qr_incld, rhofacr, nr_incld, mu_r, lamr,
        V_qr, V_nr, p3::Functions<Real,DefaultDevice>::P3Runtime());

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      crfv_device(vs).nr_incld = nr_incld[s];
      crfv_device(vs).mu_r     = mu_r[s];
      crfv_device(vs).lamr     = lamr[s];
      crfv_device(vs).V_qr     = V_qr[s];
      crfv_device(vs).V_nr     = V_nr[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(crfv_host, crfv_device);

  // Validate results
  if (SCREAM_BFB_TESTING) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(crfv_fortran[s].nr_incld == crfv_host(s).nr_incld);
      REQUIRE(crfv_fortran[s].mu_r     == crfv_host(s).mu_r);
      REQUIRE(crfv_fortran[s].lamr     == crfv_host(s).lamr);
      REQUIRE(crfv_fortran[s].V_qr     == crfv_host(s).V_qr);
      REQUIRE(crfv_fortran[s].V_nr     == crfv_host(s).V_nr);
    }
  }
}

static void run_bfb_rain_sed()
{
  auto engine = setup_random_test();

  // F90 is quite slow on weaver, so we decrease dt to reduce
  // the number of steps in rain_sed.
#ifdef EAMXX_ENABLE_GPU
  constexpr Scalar dt = 5.800E+01;
#else
  constexpr Scalar dt = 1.800E+03;
#endif

  RainSedData rsds_fortran[] = {
    //        kts, kte, ktop, kbot, kdir, dt, inv_dt, precip_liq_surf
    RainSedData(1,  72,   27,   72,   -1, dt,   1/dt,            0.0),
    RainSedData(1,  72,   72,   27,    1, dt,   1/dt,            1.0),
    RainSedData(1,  72,   27,   27,   -1, dt,   1/dt,            0.0),
    RainSedData(1,  72,   27,   27,    1, dt,   1/dt,            2.0),
  };

  static constexpr Int num_runs = sizeof(rsds_fortran) / sizeof(RainSedData);

  // Set up random input data
  for (auto& d : rsds_fortran) {
    d.randomize(engine, { {d.qr_incld, {C::QSMALL/2, C::QSMALL*2}} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  RainSedData rsds_cxx[num_runs] = {
    RainSedData(rsds_fortran[0]),
    RainSedData(rsds_fortran[1]),
    RainSedData(rsds_fortran[2]),
    RainSedData(rsds_fortran[3]),
  };

  // Get data from fortran
  for (auto& d : rsds_fortran) {
    rain_sedimentation(d);
  }

  // Get data from cxx
  for (auto& d : rsds_cxx) {
    // The code below is to force a result difference. This is used by the
    // scream/scripts internal testing to verify that various DIFFs are detected.
    auto inv_dt = d.inv_dt;
#if defined(SCREAM_FORCE_RUN_DIFF)
    inv_dt *= 2;
#endif
    rain_sedimentation_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                         d.qr_incld, d.rho, d.inv_rho, d.rhofacr, d.cld_frac_r, d.inv_dz,
                         d.dt, inv_dt,
                         d.qr, d.nr, d.nr_incld, d.mu_r, d.lamr, &d.precip_liq_surf, d.precip_liq_flux,
                         d.qr_tend, d.nr_tend);
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      // Due to pack issues, we must restrict checks to the active k space
      Int start = std::min(rsds_fortran[i].kbot, rsds_fortran[i].ktop) - 1; // 0-based indx
      Int end   = std::max(rsds_fortran[i].kbot, rsds_fortran[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(rsds_fortran[i].qr[k]              == rsds_cxx[i].qr[k]);
        REQUIRE(rsds_fortran[i].nr[k]              == rsds_cxx[i].nr[k]);
        REQUIRE(rsds_fortran[i].nr_incld[k]        == rsds_cxx[i].nr_incld[k]);
        REQUIRE(rsds_fortran[i].mu_r[k]            == rsds_cxx[i].mu_r[k]);
        REQUIRE(rsds_fortran[i].lamr[k]            == rsds_cxx[i].lamr[k]);
        REQUIRE(rsds_fortran[i].precip_liq_flux[k] == rsds_cxx[i].precip_liq_flux[k]);
        REQUIRE(rsds_fortran[i].qr_tend[k]         == rsds_cxx[i].qr_tend[k]);
        REQUIRE(rsds_fortran[i].nr_tend[k]         == rsds_cxx[i].nr_tend[k]);
      }
      REQUIRE(rsds_fortran[i].precip_liq_flux[end] == rsds_cxx[i].precip_liq_flux[end]);
      REQUIRE(rsds_fortran[i].precip_liq_surf      == rsds_cxx[i].precip_liq_surf);
    }
  }
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

  scream::p3::p3_init(); // need fortran table data

  TRS::run_phys();
  TRS::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
