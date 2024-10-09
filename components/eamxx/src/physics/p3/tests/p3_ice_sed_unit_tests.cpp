#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
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
struct UnitWrap::UnitTest<D>::TestIceSed {

static void run_phys_calc_bulk_rhime()
{
  // TODO
}

static void run_phys_ice_sed()
{
  // TODO
}

static void run_phys_homogeneous_freezing()
{
  // TODO
}

static void run_phys()
{
  run_phys_calc_bulk_rhime();
  run_phys_ice_sed();
  run_phys_homogeneous_freezing();
}

static void run_bfb_calc_bulk_rhime()
{
  constexpr Scalar qsmall = C::QSMALL;

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
  view_1d<CalcBulkRhoRimeData> cbrr_device("cbrr", max_pack_size);
  const auto cbrr_host = Kokkos::create_mirror_view(cbrr_device);
  std::copy(&cbrr_fortran[0], &cbrr_fortran[0] + max_pack_size, cbrr_host.data());
  Kokkos::deep_copy(cbrr_device, cbrr_host);

  // Get data from fortran
  for (Int i = 0; i < max_pack_size; ++i) {
    if (cbrr_fortran[i].qi_tot > qsmall) {
      calc_bulk_rho_rime(cbrr_fortran[i]);
    }
  }

  // Calc bulk rime from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack qi_tot, qi_rim, bi_rim;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      qi_tot[s] = cbrr_device(vs).qi_tot;
      qi_rim[s] = cbrr_device(vs).qi_rim;
      bi_rim[s] = cbrr_device(vs).bi_rim;
    }

    Smask gt_small(qi_tot > qsmall);
    Spack rho_rime = Functions::calc_bulk_rho_rime(
        qi_tot, qi_rim, bi_rim, p3::Functions<Real, DefaultDevice>::P3Runtime(),
        gt_small);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      cbrr_device(vs).qi_rim   = qi_rim[s];
      cbrr_device(vs).bi_rim   = bi_rim[s];
      cbrr_device(vs).rho_rime = rho_rime[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(cbrr_host, cbrr_device);

  // Validate results
  if (SCREAM_BFB_TESTING) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(cbrr_fortran[s].qi_rim   == cbrr_host(s).qi_rim);
      REQUIRE(cbrr_fortran[s].bi_rim   == cbrr_host(s).bi_rim);
      REQUIRE(cbrr_fortran[s].rho_rime == cbrr_host(s).rho_rime);
    }
  }
}

static void run_bfb_ice_sed()
{
  auto engine = setup_random_test();

  IceSedData isds_fortran[] = {
    //       kts, kte, ktop, kbot, kdir,        dt,   inv_dt, precip_ice_surf
    IceSedData(1,  72,   27,   72,   -1, 1.800E+03, 5.556E-04,            0.0),
    IceSedData(1,  72,   72,   27,    1, 1.800E+03, 5.556E-04,            1.0),
    IceSedData(1,  72,   27,   27,   -1, 1.800E+03, 5.556E-04,            0.0),
    IceSedData(1,  72,   27,   27,    1, 1.800E+03, 5.556E-04,            2.0),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(IceSedData);

  // Set up random input data
  for (auto& d : isds_fortran) {
    d.randomize(engine, { {d.qi_incld, {C::QSMALL/2, C::QSMALL*2}} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  IceSedData isds_cxx[num_runs] = {
    IceSedData(isds_fortran[0]),
    IceSedData(isds_fortran[1]),
    IceSedData(isds_fortran[2]),
    IceSedData(isds_fortran[3]),
  };

  // Get data from fortran
  for (auto& d : isds_fortran) {
    ice_sedimentation(d);
  }

  // Get data from cxx
  for (auto& d : isds_cxx) {
    ice_sedimentation_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                        d.rho, d.inv_rho, d.rhofaci, d.cld_frac_i, d.inv_dz,
                        d.dt, d.inv_dt,
                        d.qi, d.qi_incld, d.ni, d.qm, d.qm_incld, d.bm, d.bm_incld,
                        d.ni_incld, &d.precip_ice_surf, d.qi_tend, d.ni_tend);
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      // Due to pack issues, we must restrict checks to the active k space
      Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
      Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(isds_fortran[i].qi[k]       == isds_cxx[i].qi[k]);
        REQUIRE(isds_fortran[i].qi_incld[k] == isds_cxx[i].qi_incld[k]);
        REQUIRE(isds_fortran[i].ni[k]       == isds_cxx[i].ni[k]);
        REQUIRE(isds_fortran[i].ni_incld[k] == isds_cxx[i].ni_incld[k]);
        REQUIRE(isds_fortran[i].qm[k]       == isds_cxx[i].qm[k]);
        REQUIRE(isds_fortran[i].qm_incld[k] == isds_cxx[i].qm_incld[k]);
        REQUIRE(isds_fortran[i].bm[k]       == isds_cxx[i].bm[k]);
        REQUIRE(isds_fortran[i].bm_incld[k] == isds_cxx[i].bm_incld[k]);
        REQUIRE(isds_fortran[i].qi_tend[k]  == isds_cxx[i].qi_tend[k]);
        REQUIRE(isds_fortran[i].ni_tend[k]  == isds_cxx[i].ni_tend[k]);
      }
      REQUIRE(isds_fortran[i].precip_ice_surf == isds_cxx[i].precip_ice_surf);
    }
  }
}

static void run_bfb_homogeneous_freezing()
{
  constexpr Scalar latice = C::LatIce;

  auto engine = setup_random_test();

  HomogeneousFreezingData hfds_fortran[] = {
    //                    kts, kte, ktop, kbot, kdir
    HomogeneousFreezingData(1,  72,   27,   72,   -1),
    HomogeneousFreezingData(1,  72,   72,   27,    1),
    HomogeneousFreezingData(1,  72,   27,   27,   -1),
    HomogeneousFreezingData(1,  72,   27,   27,    1),
  };

  static constexpr Int num_runs = sizeof(hfds_fortran) / sizeof(HomogeneousFreezingData);

  // Set up random input data
  for (auto& d : hfds_fortran) {
    const auto qsmall_r = std::make_pair(C::QSMALL/2, C::QSMALL*2);
    d.randomize(engine, { {d.T_atm, {C::T_homogfrz - 10, C::T_homogfrz + 10}}, {d.qc, qsmall_r}, {d.qr, qsmall_r} });

    // C++ impl uses constants for latent_heat values. Manually set here
    // so F90 can match
    for (int k=0; k<d.kte; ++k) {
      d.latent_heat_fusion[k] = latice;
    }
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  HomogeneousFreezingData hfds_cxx[num_runs] = {
    HomogeneousFreezingData(hfds_fortran[0]),
    HomogeneousFreezingData(hfds_fortran[1]),
    HomogeneousFreezingData(hfds_fortran[2]),
    HomogeneousFreezingData(hfds_fortran[3]),
  };

    // Get data from fortran
  for (auto& d : hfds_fortran) {
    homogeneous_freezing(d);
  }

  // Get data from cxx
  for (auto& d : hfds_cxx) {
    homogeneous_freezing_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                           d.T_atm, d.inv_exner,
                           d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.th_atm);
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      // Due to pack issues, we must restrict checks to the active k space
      Int start = std::min(hfds_fortran[i].kbot, hfds_fortran[i].ktop) - 1; // 0-based indx
      Int end   = std::max(hfds_fortran[i].kbot, hfds_fortran[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(hfds_fortran[i].qc[k]     == hfds_cxx[i].qc[k]);
        REQUIRE(hfds_fortran[i].nc[k]     == hfds_cxx[i].nc[k]);
        REQUIRE(hfds_fortran[i].qr[k]     == hfds_cxx[i].qr[k]);
        REQUIRE(hfds_fortran[i].nr[k]     == hfds_cxx[i].nr[k]);
        REQUIRE(hfds_fortran[i].qi[k]     == hfds_cxx[i].qi[k]);
        REQUIRE(hfds_fortran[i].ni[k]     == hfds_cxx[i].ni[k]);
        REQUIRE(hfds_fortran[i].qm[k]     == hfds_cxx[i].qm[k]);
        REQUIRE(hfds_fortran[i].bm[k]     == hfds_cxx[i].bm[k]);
        REQUIRE(hfds_fortran[i].th_atm[k] == hfds_cxx[i].th_atm[k]);
      }
    }
  }
}

static void run_bfb()
{
  run_bfb_calc_bulk_rhime();
  run_bfb_ice_sed();
  run_bfb_homogeneous_freezing();
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
