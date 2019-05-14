#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include <thread>
#include <array>
#include <algorithm>

namespace {

using namespace scream;
using namespace scream::p3;

/*
 * Unit-tests for important (going to be used in scream) things.
 */

struct UnitWrap {

template <typename D=DefaultDevice>
struct UnitTest : public KokkosTypes<D> {

using Device     = D;
using MemberType = typename KokkosTypes<Device>::MemberType;
using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;
using ExeSpace   = typename KokkosTypes<Device>::ExeSpace;

template <typename S>
using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
template <typename S>
using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
template <typename S>
using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

//
// Test find_top and find_bottom
//
static void unittest_find_top_bottom()
{
  using Functions = scream::p3::Functions<Real, Device>;
  const int max_threads =
#ifdef KOKKOS_ENABLE_OPENMP
    Kokkos::OpenMP::concurrency()
#else
    1
#endif
    ;

  const int num_vert = 100;
  view_1d<Real> qr_present("qr", num_vert);
  view_1d<Real> qr_not_present("qr", num_vert);
  const int kbot = 0;
  const int ktop = num_vert - 1;
  const Real small = 0.5;
  const int large_idx_start = 33;
  const int large_idx_stop  = 77;

  auto mirror_qrp  = Kokkos::create_mirror_view(qr_present);
  auto mirror_qrnp = Kokkos::create_mirror_view(qr_not_present);

  for (int i = 0; i < num_vert; ++i) {
    mirror_qrnp(i) = small - 0.1;
    mirror_qrp(i)  = (i >= large_idx_start && i <= large_idx_stop) ? small + 0.1 : small - 0.1;
  }

  // add "hole"
  mirror_qrp(50) = small;

  Kokkos::deep_copy(qr_present, mirror_qrp);
  Kokkos::deep_copy(qr_not_present, mirror_qrnp);

  for (int team_size : {1, max_threads}) {
    const auto policy = util::ExeSpaceUtils<ExeSpace>::get_team_policy_force_team_size(1, team_size);

    int errs_for_this_ts = 0;
    Kokkos::parallel_reduce("unittest_find_top_bottom",
                            policy,
                            KOKKOS_LAMBDA(const MemberType& team, int& total_errs) {
      int nerrs_local = 0;

      //
      // Test find_top and find_bottom
      //

      bool log_qxpresent;
      int top = Functions::find_top(team, qr_present, small, kbot, ktop, 1, log_qxpresent);
      if (!log_qxpresent) ++nerrs_local;
      if (top != large_idx_stop) ++nerrs_local;

      int bot = Functions::find_bottom(team, qr_present, small, kbot, top, 1, log_qxpresent);
      if (!log_qxpresent) ++nerrs_local;
      if (bot != large_idx_start) ++nerrs_local;

      top = Functions::find_top(team, qr_present, small, ktop, kbot, -1, log_qxpresent);
      if (!log_qxpresent) ++nerrs_local;
      if (top != large_idx_start) ++nerrs_local;

      bot = Functions::find_bottom(team, qr_present, small, ktop, top, -1, log_qxpresent);
      if (!log_qxpresent) ++nerrs_local;
      if (bot != large_idx_stop) ++nerrs_local;

      top = Functions::find_top(team, qr_not_present, small, kbot, ktop, 1, log_qxpresent);
      if (log_qxpresent) ++nerrs_local;
      //if (top != 0) { std::cout << "top(" << top << ") != 0" << std::endl; ++nerrs_local; }

      bot = Functions::find_bottom(team, qr_not_present, small, kbot, ktop, 1, log_qxpresent);
      if (log_qxpresent) ++nerrs_local;
      //if (bot != 0) { std::cout << "bot(" << bot << ") != 0" << std::endl; ++nerrs_local; }

      total_errs += nerrs_local;
    }, errs_for_this_ts);

    REQUIRE(errs_for_this_ts == 0);
  }
}

// Derive rough bounds on table domain:
//
// First, mu_r is in [0,9].
//
// Second, start with
//     dum1 = (1 + mu_r)/lamr
// Let's get the lower limit on lamr as a function of mu first. So
//     1 = rdumii = (dum1 1e6 + 5) inv_dum3, inv_dum3 = 0.1
//       = ((mu_r+1)/lamr 1e6 + 5) 0.1
//     (1/0.1 - 5) 1e-6 lamr = 1 + mu_r
//     lamr = 0.2e6 (1 + mu_r).
// Now the upper limit:
//     300 = (dum1 1e6 + 5) inv_dum3, inv_dum3 = 1/30
//         = ((mu_r+1)/lamr 1e6 - 195) 1/30 + 20
//     => lamr = 1e6 (1 + mu_r) / (30 (300 - 20) + 195)
//             ~ 116.3 (1 + mu_r)
//
// This test checks for smoothness (C^0 function) of the interpolant both inside
// the domains over which the tables are defined and outside them. The primary
// tool is measuring the maximum slope magnitude as a function of mesh
// refinement, where the mesh is a 1D mesh transecting the table domain.

struct TestTable3 {
  using Functions = scream::p3::Functions<Real, Device>;
  using view_1d_table = typename Functions::view_1d_table;
  using view_2d_table = typename Functions::view_2d_table;
  using Scalar = typename Functions::Scalar;
  using Smask = typename Functions::Smask;
  using Spack = typename Functions::Spack;
  using Table3 = typename Functions::Table3;
  using RangePolicy = Kokkos::RangePolicy<typename Device::execution_space>;

  KOKKOS_FUNCTION static Scalar calc_lamr (const Scalar& mu_r, const Scalar& alpha) {
    // Parameters for lower and upper bounds, derived above, multiplied by
    // factors so we go outside of the bounds a bit. Using these, we map alpha
    // in [0,1] -> meaningful lamr.
    const Scalar lamr_lo = 0.1*116.3, lamr_hi = 1.5*0.2e6;
    return (lamr_lo + alpha*(lamr_hi - lamr_lo))*(1 + mu_r);
  }

  KOKKOS_FUNCTION static Scalar calc_mu_r (const Scalar& alpha) {
    return alpha*10;
  }

  // Perform the table lookup and interpolation operations for (mu_r, lamr).
  KOKKOS_FUNCTION static Spack interp (const view_2d_table& table, const Scalar& mu_r,
                                       const Scalar& lamr) {
    // Init the pack to all the same value, and compute in every pack slot.
    Smask qr_gt_small(true);
    Spack mu_r_p(mu_r), lamr_p(lamr);
    Table3 t3;
    Functions::lookup(qr_gt_small, mu_r_p, lamr_p, t3);
    Spack val(qr_gt_small, Functions::apply_table(qr_gt_small, table, t3));
    return val;
  }

  static void run () {
    // This test doesn't use mu_r_table, as that is not a table3 type. It
    // doesn't matter whether we use vm_table or vn_table, as the table values
    // don't matter in what we are testing; we are testing interpolation
    // procedures that are indepennt of particular table values.
    view_1d_table mu_r_table;
    view_2d_table vn_table, vm_table;
    Functions::init_kokkos_tables(vn_table, vm_table, mu_r_table);

    // Estimate two maximum slope magnitudes for two meshes, the second 10x
    // refined w.r.t. the first.
    Real slopes[2];
    const Int nslopes = sizeof(slopes)/sizeof(*slopes);
    Int N;

    // Study the mu_r direction.
    // For a sequence of refined meshes:
    N = 1000;
    for (Int refine = 0; refine < nslopes; ++refine) {
      // Number of cells in the mesh.
      N *= 10;
      // Cell size relative to a parameter domain of 1.
      const Scalar delta = 1.0/N;

      // Compute the slope magnitude at a specific (mu_r, lamr) in the mu_r
      // direction.
      const Scalar lamr = calc_lamr(4.5, 0.5);
      const auto get_max_slope = KOKKOS_LAMBDA (const Int& i, Scalar& slope) {
        // Interpolate at a specific (mu_r, lamr).
        const auto eval = [&] (const Int& i) {
          const auto alpha = double(i)/N;
          const auto mu_r = calc_mu_r(alpha);
          const auto val = interp(vm_table, mu_r, lamr);
          return std::log(val[0]);
        };
        slope = util::max(slope, std::abs((eval(i+1) - eval(i))/delta));
      };

      Scalar max_slope;
      Kokkos::parallel_reduce(RangePolicy(0, N), get_max_slope,
                              Kokkos::Max<Scalar>(max_slope));
      Kokkos::fence();
      slopes[refine] = max_slope;
    }

    // Now that we have collected slopes as a function of 10x mesh refinement,
    // determine whether the slope estimates are converging. If they are, we can
    // conclude there are no discontinuities. If they are not, then we are sensing
    // a discontinuity, which is a bug.
    //   In detail, for a 10x mesh refinement, a good slope growth rate is right
    // around 1, and a bad one is is roughly 10. We set the threshold at 1.1.
    const auto check_growth = [&] (const std::string& label, const Scalar& growth) {
      if (growth > 1.1) {
        std::cout << "Table3 FAIL: Slopes in the " << label << " direction are "
        << slopes[0] << " and " << slopes[1]
        << ", which grows by factor " << growth
        << ". Near 1 is good; near 10 is bad.\n";
        REQUIRE(false);
      }
    };
    check_growth("mu_r", slopes[1]/slopes[0]);

    // Study the lamr direction.
    N = 4000;
    for (Int refine = 0; refine < nslopes; ++refine) {
      N *= 2;
      const Scalar delta = 1.0/N;

      // Compute the slope magnitude at a specific (mu_r, lamr) in the lamr
      // direction.
      const Scalar mu_r = 3.5;
      const auto get_max_slope = KOKKOS_LAMBDA (const Int& i, Scalar& slope) {
        // Interpolate at a specific (mu_r, lamr).
        const auto eval = [&] (const Int& i) {
          const auto alpha = double(i)/N;
          const auto lamr = calc_lamr(mu_r, alpha);
          const auto val = interp(vm_table, mu_r, lamr);
          return std::log(val[0]);
        };
        slope = util::max(slope, std::abs((eval(i+1) - eval(i))/delta));
      };

      Scalar max_slope;
      Kokkos::parallel_reduce(RangePolicy(0, N), get_max_slope,
                              Kokkos::Max<Scalar>(max_slope));
      Kokkos::fence();
      slopes[refine] = max_slope;
    }
    check_growth("lamr", slopes[1]/slopes[0]);
  }
};

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
static void unittest_upwind () {
  static const Int nfield = 2;

  using Functions = scream::p3::Functions<Real, Device>;
  using Scalar = typename Functions::Scalar;
  using Pack = typename Functions::Pack;
  using Spack = typename Functions::Spack;

  const auto eps = std::numeric_limits<Scalar>::epsilon();

  Int nerr = 0;
  for (Int nk : {17, 32, 77, 128}) {
    const Int npack = (nk + Pack::n - 1) / Pack::n, kmin = 0, kmax = nk - 1;
    const Real max_speed = 4.2, min_dz = 0.33;
    const Real dt = min_dz/max_speed;

    view_1d<Pack> rho("rho", npack), inv_rho("inv_rho", npack), inv_dz("inv_dz", npack);
    const auto lrho = smallize(rho), linv_rho = smallize(inv_rho), linv_dz = smallize(inv_dz);

    Kokkos::Array<view_1d<Pack>, nfield> flux, V, r;
    Kokkos::Array<ko::Unmanaged<view_1d<Spack> >, nfield> lflux, lV, lr;
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

};
};

TEST_CASE("p3_find", "[p3_functions]")
{
  UnitWrap::UnitTest<scream::DefaultDevice>::unittest_find_top_bottom();
}

TEST_CASE("p3_tables", "[p3_functions]")
{
  // Populate tables with contrived values, don't need realistic values
  // for unit testing
  using Globals = scream::p3::Globals<Real>;
  for (size_t i = 0; i < Globals::VN_TABLE.size(); ++i) {
    for (size_t j = 0; j < Globals::VN_TABLE[0].size(); ++j) {
      Globals::VN_TABLE[i][j] = 1 + i + j;
      Globals::VM_TABLE[i][j] = 1 + i * j;
    }
  }

  for (size_t i = 0; i < Globals::MU_R_TABLE.size(); ++i) {
    Globals::MU_R_TABLE[i] = 1 + i;
  }
  UnitWrap::UnitTest<scream::DefaultDevice>::TestTable3::run();
}

TEST_CASE("p3_upwind", "[p3_functions]")
{
  UnitWrap::UnitTest<scream::DefaultDevice>::unittest_upwind();
}

} // namespace
