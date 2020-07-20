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

// r[1] is a mixing ratio field. The test advects r[1] some number of time
// steps. The test checks global mass conservation and extrema non-violation at
// each step. B/c of the consistency issue noted in the
// calc_first_order_upwind_step doc, r[0] = 1 initially so that r[0]*rho is the
// true, i.e. correctly advected, total density. Note that r[0] will not remain
// uniformly 1. Extrema of r[1]/r[0], the true mixing ratio, at time step n+1
// must be within the min and max of r[1]/r[0] at time step n. Mass conservation
// includes the material fluxed out of the domain. r[1] must be initialized with
// an initial condition. The details of the IC are not important except that the
// profile should be nontrivial. Also, it is initialized so the first and last
// cells in the domain are 0. This lets us check the restricted-domain usage of
// the upwind routine in the first time step.
template <typename D>
struct UnitWrap::UnitTest<D>::TestUpwind {

static void run_phys()
{
  static const Int nfield = 2;

  const auto eps = std::numeric_limits<Scalar>::epsilon();

  Int nerr = 0;
  for (Int nk : {17, 32, 77, 128}) {
    const Int npack = (nk + Pack::n - 1) / Pack::n, kmin = 0, kmax = nk - 1;
    const Real max_speed = 4.2, min_dz = 0.33;
    const Real dt = min_dz/max_speed;

    view_1d<Pack> rho("rho", npack), inv_rho("inv_rho", npack), inv_dz("inv_dz", npack);
    const auto lrho = smallize(rho), linv_rho = smallize(inv_rho), linv_dz = smallize(inv_dz);

    Kokkos::Array<view_1d<Pack>, nfield> flux, V, r;
    Kokkos::Array<uview_1d<Spack>, nfield> lflux, lV, lr;
    const auto init_array = [&] (const std::string& name, const Int& i, decltype(flux)& f,
                                 decltype(lflux)& lf) {
      f[i] = view_1d<Pack>("f", npack);
      lf[i] = smallize(f[i]);
    };
    for (int i = 0; i < nfield; ++i) {
      init_array("flux", i, flux, lflux);
      init_array("V", i, V, lV);
      init_array("r", i, r, lr);
    }

    for (Int kdir : {-1, 1}) {
      const Int k_bot = kdir == 1 ? kmin : kmax;
      const Int k_top = kdir == 1 ? kmax : kmin;

      // Set rho, dz, mixing ratio r.
      const auto init_fields = KOKKOS_LAMBDA (const MemberType& team) {
        const auto set_fields = [&] (const Int& k) {
          for (Int i = 0; i < nfield; ++i) {
            const auto range = scream::pack::range<Pack>(k*Pack::n);
            rho(k) = 1 + range/nk;
            inv_rho(k) = 1 / rho(k);
            inv_dz(k) = 1 / (min_dz + range*range / (nk*nk));
            V[i](k) = 0.5*(1 + range/nk) * max_speed;
            if (i == 1) {
              r[i](k) = 0;
              const auto mask = range >= 2 && range < nk-2;
              r[i](k).set(mask, range/nk); // Nontrivial mixing ratio.
            } else {
              r[i](k) = 1; // Evolve the background density field.
            }
            scream_kassert((V[i](k) >= 0).all());
            scream_kassert((V[i](k) <= max_speed || (range >= nk)).all());
          }
          scream_kassert((V[0](k) == V[1](k)).all());
        };
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, npack), set_fields);
        team.team_barrier();
      };
      Kokkos::parallel_for(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, npack),
                           init_fields);

      const auto sflux = scalarize(flux[1]);
      for (Int time_step = 0; time_step < 2*nk; ++time_step) {
        // Take one upwind step.
        const auto step = KOKKOS_LAMBDA (const MemberType& team, Int& nerr) {
          const auto sr = scalarize(r[1]), srho = scalarize(rho), sinv_dz = scalarize(inv_dz);
          const auto sr0 = scalarize(r[0]);

          // Gather diagnostics: total mass and extremal mixing ratio values.
          const auto gather_diagnostics = [&] (Scalar& mass, Scalar& r_min, Scalar& r_max) {
            mass = 0;
            const auto sum_mass = [&] (const Int& k, Scalar& mass) {
              mass += srho(k)*sr(k)/sinv_dz(k);
            };
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nk), sum_mass, mass);

            const auto find_max_r = [&] (const Int& k, Scalar& r_max) {
              // The background rho is is not advected in P3. Thus, here we
              // advect the initially uniform mixing ratio r[0] to capture the
              // true advected background rho_true = r[0]*rho. Then the true
              // mixing ratio corresponding to r[1] is
              //     r_true = (r[1]*rho)/(rho_true)
              //            = (r[1]*rho)/(r[0]*rho) = r[1]/r[0].
              // This mixing ratio is tested to show that it does not violate
              // the previous time step's global extrema.
              //
              // At the inflow boundary, where inflow is 0, we need to be
              // careful about sr(k_top)/sr0(k_top) becoming dominated by noise
              // as each goes to 0. eps^2 relative to a starting value of
              // sr0(k_top) = 1 at time 0 is unnecessarily small (we could
              // choose a larger lower bound and still be testing things well),
              // but it works, so we might as well use it.
              if (eps*sr0(k) < eps) return;
              const auto mixing_ratio_true = sr(k)/sr0(k);
              r_max = util::max(mixing_ratio_true, r_max);
            };
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nk), find_max_r,
                                    Kokkos::Max<Scalar>(r_max));

            const auto find_min_r = [&] (const Int& k, Scalar& r_min) {
              const auto mixing_ratio_true = sr(k)/sr0(k);
              r_min = util::min(mixing_ratio_true, r_min);
            };
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nk), find_min_r,
                                    Kokkos::Min<Scalar>(r_min));
          };

          // Gather diagnostics before the step.
          Scalar mass0, r_min0, r_max0;
          gather_diagnostics(mass0, r_min0, r_max0);
          team.team_barrier();

          // Compute.
          Int k_bot_lcl = k_bot, k_top_lcl = k_top;
          if (time_step == 0) {
            // In the first step, the IC for r[1] is such that we can test
            // restricting the interval. But the IC for r[0] does not permit
            // it. Thus, make two calls to the upwind routine:
            //   1. Full domain for r[0].
            Functions::calc_first_order_upwind_step(
              lrho, linv_rho, linv_dz, team, nk, k_bot, k_top, kdir, dt,
              lflux[0], lV[0], lr[0]);
            k_bot_lcl += kdir;
            k_top_lcl -= kdir;
            //   2. Restricted domain for r[1] in first time step only. Note
            // that the restriction is unnecesary but just here to test the
            // restriction code.
            Functions::calc_first_order_upwind_step(
              lrho, linv_rho, linv_dz, team, nk, k_bot_lcl, k_top_lcl, kdir, dt,
              lflux[1], lV[1], lr[1]);
          } else {
            Functions::template calc_first_order_upwind_step<nfield>(
              lrho, linv_rho, linv_dz, team, nk, k_bot_lcl, k_top_lcl, kdir, dt,
              {&lflux[0], &lflux[1]}, {&lV[0], &lV[1]}, {&lr[0], &lr[1]});
          }
          team.team_barrier();

          // Gather diagnostics after the step.
          Scalar mass1, r_min1, r_max1;
          gather_diagnostics(mass1, r_min1, r_max1);
          // Include mass flowing out of the boundary.
          if (time_step > 1) mass1 += sflux(kdir == 1 ? 0 : nk-1)*dt;
          team.team_barrier();

          // Check diagnostics.
          //   1. Check for conservation of mass.
          if (util::reldif(mass0, mass1) > 1e1*eps) ++nerr;
          //   2. Check for non-violation of global extrema.
          if (r_min1 < r_min0 - 10*eps) ++nerr;
          if (r_max1 > r_max0 + 10*eps) ++nerr;
        };
        Int lnerr;
        Kokkos::parallel_reduce(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, npack),
                                step, lnerr);
        nerr += lnerr;
        Kokkos::fence();
        REQUIRE(nerr == 0);
      }
    }
  }
}

static void run_bfb()
{
  CalcUpwindData cuds_fortran[] = {
                // kts, kte, kdir, kbot, k_qxtop, na,   dt_sub,  rho range, inv_dzq range, vs range, qnx range
    CalcUpwindData(  1,  72,   -1,   72,      36,  2,  1.833E+03,
                   std::make_pair(4.056E-03, 1.153E+00),
                   std::make_pair(2.863E-05, 8.141E-03),
                   std::make_pair(2.965E-02, 3.555E+00),
                   std::make_pair(7.701E-16, 2.119E-04)),

    CalcUpwindData(  1,  72,    1,   36,      72,  2,  1.833E+03,
                   std::make_pair(4.056E-03, 1.153E+00),
                   std::make_pair(2.863E-05, 8.141E-03),
                   std::make_pair(2.965E-02, 3.555E+00),
                   std::make_pair(7.701E-16, 2.119E-04)),

    CalcUpwindData(  1,  72,   -1,   72,      36,  4,  1.833E+03,
                   std::make_pair(4.056E-03, 1.153E+00),
                   std::make_pair(2.863E-05, 8.141E-03),
                   std::make_pair(2.965E-02, 3.555E+00),
                   std::make_pair(7.701E-16, 2.119E-04)),

    CalcUpwindData(  1,  72,   -1,   72,      72,  2,  1.833E+03,
                   std::make_pair(4.056E-03, 1.153E+00),
                   std::make_pair(2.863E-05, 8.141E-03),
                   std::make_pair(2.965E-02, 3.555E+00),
                   std::make_pair(7.701E-16, 2.119E-04)),

    CalcUpwindData(  1,  32,   -1,   24,      8,  2,  1.833E+03,
                   std::make_pair(4.056E-03, 1.153E+00),
                   std::make_pair(2.863E-05, 8.141E-03),
                   std::make_pair(2.965E-02, 3.555E+00),
                   std::make_pair(7.701E-16, 2.119E-04)),

    CalcUpwindData(  1,  32,    1,    7,      21,  1,  1.833E+03,
                   std::make_pair(4.056E-03, 1.153E+00),
                   std::make_pair(2.863E-05, 8.141E-03),
                   std::make_pair(2.965E-02, 3.555E+00),
                   std::make_pair(7.701E-16, 2.119E-04)),

    CalcUpwindData(  1,  32,   -1,   21,      7,  1,  1.833E+03,
                   std::make_pair(4.056E-03, 1.153E+00),
                   std::make_pair(2.863E-05, 8.141E-03),
                   std::make_pair(2.965E-02, 3.555E+00),
                   std::make_pair(7.701E-16, 2.119E-04)),

  };

  static constexpr Int num_runs = sizeof(cuds_fortran) / sizeof(CalcUpwindData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  CalcUpwindData cuds_cxx[num_runs] = {
    CalcUpwindData(cuds_fortran[0]),
    CalcUpwindData(cuds_fortran[1]),
    CalcUpwindData(cuds_fortran[2]),
    CalcUpwindData(cuds_fortran[3]),
    CalcUpwindData(cuds_fortran[4]),
    CalcUpwindData(cuds_fortran[5]),
    CalcUpwindData(cuds_fortran[6]),
  };

  // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    calc_first_order_upwind_step(cuds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    calc_first_order_upwind_step_f(
      cuds_cxx[i].kts, cuds_cxx[i].kte, cuds_cxx[i].kdir, cuds_cxx[i].kbot, cuds_cxx[i].k_qxtop, cuds_cxx[i].dt_sub,
      cuds_cxx[i].rho, cuds_cxx[i].inv_rho, cuds_cxx[i].inv_dzq,
      cuds_cxx[i].num_arrays, cuds_cxx[i].fluxes, cuds_cxx[i].vs, cuds_cxx[i].qnx);
  }

  for (Int i = 0; i < num_runs; ++i) {
    for (int n = 0; n < cuds_fortran[i].num_arrays; ++n) {
      // Due to pack issues, we must restrict checks to the active k space
      Int start = std::min(cuds_fortran[i].kbot, cuds_fortran[i].k_qxtop) - 1; // 0-based indx
      Int end   = std::max(cuds_fortran[i].kbot, cuds_fortran[i].k_qxtop); // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(cuds_fortran[i].fluxes[n][k] == cuds_cxx[i].fluxes[n][k]);
        REQUIRE(cuds_fortran[i].qnx[n][k]    == cuds_cxx[i].qnx[n][k]);
      }
    }
  }
}

};

template <typename D>
struct UnitWrap::UnitTest<D>::TestGenSed {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
//   Co_max = 9.196837456784E-02 - 9.827330928362E+01
//     dt_left =1.818181818182E+01 - 1.800000000000E+03
// prt_accum = 4.959754212038E-05 - 6.211745579368E-07

  // GenSedData(Int kts_, Int kte_, Int kdir_, Int k_qxtop_, Int k_qxbot_, Int kbot_, Real Co_max_, Real dt_left_,
  //            Real prt_accum_, Int num_arrays_,
  //            std::pair<Real, Real> rho_range, std::pair<Real, Real> inv_dzq_range,
  //            std::pair<Real, Real> vs_range, std::pair<Real, Real> qnx_range);

  GenSedData gsds_fortran[] = {
    //       kts, kte, kdir, k_qxtop, k_qxbot, kbot,     Co_max,   dt_left, prt_accum, num_arrays
    GenSedData(1,  72,    -1,     36,      72,   72,  9.196E-02, 1.818E+01, 4.959E-05, 2,
               std::make_pair(4.056E-03, 1.153E+00),
               std::make_pair(2.863E-05, 8.141E-03),
               std::make_pair(2.965E-02, 3.555E+00),
               std::make_pair(7.701E-16, 2.119E-04)),

    GenSedData(1,  72,    -1,     36,      57,   72,  4.196E-01, 1.418E+02, 4.959E-06, 1,
               std::make_pair(4.056E-03, 1.153E+00),
               std::make_pair(2.863E-05, 8.141E-03),
               std::make_pair(2.965E-02, 3.555E+00),
               std::make_pair(7.701E-16, 2.119E-04)),

    GenSedData(1,  72,     1,     57,      37,   36,  4.196E-01, 1.418E+02, 4.959E-06, 1,
               std::make_pair(4.056E-03, 1.153E+00),
               std::make_pair(2.863E-05, 8.141E-03),
               std::make_pair(2.965E-02, 3.555E+00),
               std::make_pair(7.701E-16, 2.119E-04)),

    GenSedData(1,  72,    -1,     72,      72,   72,  4.196E-01, 1.418E+02, 4.959E-06, 1,
               std::make_pair(4.056E-03, 1.153E+00),
               std::make_pair(2.863E-05, 8.141E-03),
               std::make_pair(2.965E-02, 3.555E+00),
               std::make_pair(7.701E-16, 2.119E-04)),

  };

  static constexpr Int num_runs = sizeof(gsds_fortran) / sizeof(GenSedData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  GenSedData gsds_cxx[num_runs] = {
    GenSedData(gsds_fortran[0]),
    GenSedData(gsds_fortran[1]),
    GenSedData(gsds_fortran[2]),
    GenSedData(gsds_fortran[3]),
  };

  // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    generalized_sedimentation(gsds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    generalized_sedimentation_f(gsds_cxx[i].kts, gsds_cxx[i].kte, gsds_cxx[i].kdir, gsds_cxx[i].k_qxtop,
                                &gsds_cxx[i].k_qxbot, gsds_cxx[i].kbot, gsds_cxx[i].Co_max,
                                &gsds_cxx[i].dt_left, &gsds_cxx[i].prt_accum, gsds_cxx[i].inv_dzq, gsds_cxx[i].inv_rho, gsds_cxx[i].rho,
                                gsds_cxx[i].num_arrays, gsds_cxx[i].vs, gsds_cxx[i].fluxes, gsds_cxx[i].qnx);
  }

  for (Int i = 0; i < num_runs; ++i) {
    for (int n = 0; n < gsds_fortran[i].num_arrays; ++n) {
      // Due to pack issues, we must restrict checks to the active k space
      Int start = std::min(gsds_fortran[i].k_qxbot, gsds_fortran[i].k_qxtop) - 1; // 0-based indx
      Int end   = std::max(gsds_fortran[i].k_qxbot, gsds_fortran[i].k_qxtop); // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(gsds_fortran[i].fluxes[n][k] == gsds_cxx[i].fluxes[n][k]);
        REQUIRE(gsds_fortran[i].qnx[n][k]    == gsds_cxx[i].qnx[n][k]);
      }
      REQUIRE(gsds_fortran[i].k_qxbot   == gsds_cxx[i].k_qxbot);
      REQUIRE(gsds_fortran[i].dt_left   == gsds_cxx[i].dt_left);
      REQUIRE(gsds_fortran[i].prt_accum == gsds_cxx[i].prt_accum);
    }
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_upwind", "[p3_functions]")
{
  using TU = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUpwind;

  TU::run_phys();
  TU::run_bfb();
}

TEST_CASE("p3_gen_sed", "[p3_functions]")
{
  using TG = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGenSed;

  TG::run_phys();
  TG::run_bfb();
}

} // namespace
