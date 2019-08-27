#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace {

using namespace scream;
using namespace scream::p3;

/*
 * Unit-tests for p3_functions.
 */

struct UnitWrap {

template <typename D=DefaultDevice>
struct UnitTest : public KokkosTypes<D> {

using Device      = D;
using MemberType  = typename KokkosTypes<Device>::MemberType;
using TeamPolicy  = typename KokkosTypes<Device>::TeamPolicy;
using RangePolicy = typename KokkosTypes<Device>::RangePolicy;
using ExeSpace    = typename KokkosTypes<Device>::ExeSpace;

template <typename S>
using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
template <typename S>
using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
template <typename S>
using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

using Functions          = scream::p3::Functions<Real, Device>;
using view_itab_table    = typename Functions::view_itab_table;
using view_itabcol_table = typename Functions::view_itabcol_table;
using view_1d_table      = typename Functions::view_1d_table;
using view_2d_table      = typename Functions::view_2d_table;
using Scalar             = typename Functions::Scalar;
using Spack              = typename Functions::Spack;
using Pack               = typename Functions::Pack;
using IntSmallPack       = typename Functions::IntSmallPack;
using Smask              = typename Functions::Smask;
using TableIce           = typename Functions::TableIce;
using TableRain          = typename Functions::TableRain;
using Table3             = typename Functions::Table3;
using C                  = typename Functions::C;

//
// Test find_top and find_bottom
//
static void unittest_find_top_bottom()
{
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

struct TestTableIce {

  static void test_read_lookup_tables_bfb()
  {
    // Read in ice tables
    view_itab_table itab;
    view_itabcol_table itabcol;
    Functions::init_kokkos_ice_lookup_tables(itab, itabcol);

    // Get data from fortran
    P3InitAFortranData d;
    p3_init_a(d);

    // Copy device data to host
    const auto itab_host = Kokkos::create_mirror_view(itab);
    const auto itabcol_host = Kokkos::create_mirror_view(itabcol);
    Kokkos::deep_copy(itab_host, itab);
    Kokkos::deep_copy(itabcol_host, itabcol);

    // Compare (on host)
    for (size_t i = 0; i < itab_host.extent(0); ++i) {
      for (size_t j = 0; j < itab_host.extent(1); ++j) {
        for (size_t k = 0; k < itab_host.extent(2); ++k) {

          for (size_t l = 0; l < itab_host.extent(3); ++l) {
            REQUIRE(itab_host(i, j, k, l) == d.itab(i, j, k, l));
          }

          for (size_t l = 0; l < itabcol_host.extent(3); ++l) {
            for (size_t m = 0; m < itabcol_host.extent(4); ++m) {
              REQUIRE(itabcol_host(i, j, k, l, m) == d.itabcol(i, j, k, l, m));
            }
          }

        }
      }
    }
  }

  template <typename View>
  static void init_table_linear_dimension(View& table, int linear_dimension)
  {
    // set up views
    using NonConstView = typename View::non_const_type;
    const auto view_device = NonConstView("non const view");
    const auto view_host   = Kokkos::create_mirror_view(view_device);

    std::default_random_engine generator;
    std::uniform_real_distribution<Real> val_dist(0.0,100.0);

    // populate lin-dim-0 with random values, make sure values are linear
    // in the linear_dimension
    for (size_t i = 0; i < table.extent(0); ++i) {
      for(size_t j = 0; j < table.extent(1); ++j) {
        for (size_t k = 0; k < table.extent(2); ++k) {
          for (size_t l = 0; l < table.extent(3); ++l) {
            size_t dims[] = {i, j, k, l};
            if (dims[linear_dimension] == 0) {
              view_host(i, j, k, l) = val_dist(generator);
            }
            else {
              dims[linear_dimension] -= 1;
              view_host(i, j, k, l) = view_host(dims[0], dims[1], dims[2], dims[3]) + 1.0;
            }
          }
        }
      }
    }

    // Copy back to device
    Kokkos::deep_copy(view_device, view_host);
    table = view_device;
  }

  static void run_bfb()
  {
    // Read in ice tables
    view_itab_table itab;
    view_itabcol_table itabcol;
    Functions::init_kokkos_ice_lookup_tables(itab, itabcol);

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    // Load some lookup inputs, need at least one per pack value
    LookupIceData lid[max_pack_size] = {
      {0.971E-07, 0.657E+06, 0.971E-07, 0.900E+03},
      {0.510E-02, 0.454E+06, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+06, 0.811E-10, 0.500E+02},

      {0.971E-07, 0.657E+06, 0.971E-07, 0.900E+03},
      {0.510E-02, 0.454E+06, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+06, 0.811E-10, 0.500E+02},

      {0.971E-07, 0.657E+06, 0.971E-07, 0.900E+03},
      {0.510E-02, 0.454E+06, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+06, 0.811E-10, 0.500E+02},

      {0.971E-07, 0.657E+06, 0.971E-07, 0.900E+03},
      {0.510E-02, 0.454E+06, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+06, 0.811E-10, 0.500E+02}
    };

    LookupIceDataB lidb[max_pack_size] = {
      {0.263E-05, 0.100E+07},
      {0.100E-01, 0.100E+07},
      {0.000E+00, 0.100E-15},
      {0.263E-05, 0.100E+07},

      {0.263E-05, 0.100E+07},
      {0.100E-01, 0.100E+07},
      {0.000E+00, 0.100E-15},
      {0.263E-05, 0.100E+07},

      {0.263E-05, 0.100E+07},
      {0.100E-01, 0.100E+07},
      {0.000E+00, 0.100E-15},
      {0.263E-05, 0.100E+07},

      {0.263E-05, 0.100E+07},
      {0.100E-01, 0.100E+07},
      {0.000E+00, 0.100E-15},
      {0.263E-05, 0.100E+07}
    };

    static constexpr Int access_table_index = 2;
    AccessLookupTableData altd[max_pack_size] = {
      {lid[0], access_table_index},
      {lid[1], access_table_index},
      {lid[2], access_table_index},
      {lid[3], access_table_index},

      {lid[4], access_table_index},
      {lid[5], access_table_index},
      {lid[6], access_table_index},
      {lid[7], access_table_index},

      {lid[8], access_table_index},
      {lid[9], access_table_index},
      {lid[10], access_table_index},
      {lid[11], access_table_index},

      {lid[12], access_table_index},
      {lid[13], access_table_index},
      {lid[14], access_table_index},
      {lid[15], access_table_index}
    };

    AccessLookupTableCollData altcd[max_pack_size] = {
      {lid[0], lidb[0], access_table_index},
      {lid[1], lidb[1], access_table_index},
      {lid[2], lidb[2], access_table_index},
      {lid[3], lidb[3], access_table_index},

      {lid[4], lidb[4], access_table_index},
      {lid[5], lidb[5], access_table_index},
      {lid[6], lidb[6], access_table_index},
      {lid[7], lidb[7], access_table_index},

      {lid[8], lidb[8], access_table_index},
      {lid[9], lidb[9], access_table_index},
      {lid[10], lidb[10], access_table_index},
      {lid[11], lidb[11], access_table_index},

      {lid[12], lidb[12], access_table_index},
      {lid[13], lidb[13], access_table_index},
      {lid[14], lidb[14], access_table_index},
      {lid[15], lidb[15], access_table_index}
    };

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      find_lookuptable_indices_1a(lid[i]);
      find_lookuptable_indices_1b(lidb[i]);
      access_lookup_table(altd[i]);
      access_lookup_table_coll(altcd[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    view_1d<IntSmallPack> int_results("int results", 5);
    view_1d<Spack> real_results("real results", 7);
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      Smask qiti_gt_small(true);

      // Init packs
      TableIce ti;
      TableRain tr;
      Spack qitot, nitot, qirim, rhop, qr, nr;
      for(Int s = 0; s < Spack::n; ++s) {
        qitot[s] = lid[s].qitot;
        nitot[s] = lid[s].nitot;
        qirim[s] = lid[s].qirim;
        rhop[s]  = lid[s].rhop;

        qr[s]    = lidb[s].qr;
        nr[s]    = lidb[s].nr;
      }

      Functions::lookup_ice(qiti_gt_small, qitot, nitot, qirim, rhop, ti);
      Functions::lookup_rain(qiti_gt_small, qr, nr, tr);
      Spack ice_result = Functions::apply_table_ice(qiti_gt_small, access_table_index-1, itab, ti);
      Spack rain_result = Functions::apply_table_coll(qiti_gt_small, access_table_index-1, itabcol, ti, tr);

      int_results(0) = ti.dumi;
      int_results(1) = ti.dumjj;
      int_results(2) = ti.dumii;
      int_results(3) = ti.dumzz;

      int_results(4) = tr.dumj;

      real_results(0) = ti.dum1;
      real_results(1) = ti.dum4;
      real_results(2) = ti.dum5;
      real_results(3) = ti.dum6;

      real_results(4) = tr.dum3;

      real_results(5) = ice_result;

      real_results(6) = rain_result;
    });
    auto int_results_mirror  = Kokkos::create_mirror_view(int_results);
    auto real_results_mirror = Kokkos::create_mirror_view(real_results);
    Kokkos::deep_copy(int_results_mirror, int_results);
    Kokkos::deep_copy(real_results_mirror, real_results);

    // Validate results
    for(int s = 0; s < Spack::n; ++s) {
      // +1 for O vs 1-based indexing
      REQUIRE(int_results_mirror(0)[s]+1 == lid[s].dumi);
      REQUIRE(int_results_mirror(1)[s]+1 == lid[s].dumjj);
      REQUIRE(int_results_mirror(2)[s]+1 == lid[s].dumii);
      REQUIRE(int_results_mirror(3)[s]+1 == lid[s].dumzz);

      REQUIRE(int_results_mirror(4)[s]+1 == lidb[s].dumj);

      REQUIRE(real_results_mirror(0)[s] == lid[s].dum1);
      REQUIRE(real_results_mirror(1)[s] == lid[s].dum4);
      REQUIRE(real_results_mirror(2)[s] == lid[s].dum5);
      REQUIRE(real_results_mirror(3)[s] == lid[s].dum6);

      REQUIRE(real_results_mirror(4)[s] == lidb[s].dum3);

      REQUIRE(real_results_mirror(5)[s] == altd[s].proc);

      REQUIRE(real_results_mirror(6)[s] == altcd[s].proc);
    }
  }

  static void run_phys()
  {
    view_itab_table itab;
    init_table_linear_dimension(itab, 0);

    int nerr = 0;
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(itab.extent(0), itab.extent(1)));
    Kokkos::parallel_reduce("TestTableIce::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {
      //int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, itab.extent(1)), [&] (const int& j) {

        for (size_t k = 0; k < itab.extent(2); ++k) {
          for (size_t l = 0; l < itab.extent(3); ++l) {
            Smask qiti_gt_small(true);

            // Init packs to same value, TODO: how to pick use values?
            Spack qitot(0.1), nitot(0.2), qirim(0.3), rhop(0.4), qr(0.5), nr(0.6);

            TableIce ti;
            TableRain tr;
            Functions::lookup_ice(qiti_gt_small, qitot, nitot, qirim, rhop, ti);
            Functions::lookup_rain(qiti_gt_small, qr, nr, tr);

            /*Spack proc1 = */ Functions::apply_table_ice(qiti_gt_small, 1, itab, ti);
            //Spack proc2 = Functions::apply_table_coll(qiti_gt_small, 1, itabcol, ti, tr);

            // TODO: how to test?
          }
        }
      });
      errors = 0;

    }, nerr);

    view_itabcol_table itabcol;

    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
};

struct TestP3Func
{

  KOKKOS_FUNCTION  static void saturation_tests(const Scalar& temperature, const Scalar& pressure, const Scalar& correct_sat_ice_p,
    const Scalar& correct_sat_liq_p, const Scalar&  correct_mix_ice_r, const Scalar& correct_mix_liq_r, int& errors ){

    const Spack temps(temperature);
    const Spack pres(pressure);

    Spack sat_ice_p = Functions::polysvp1(temps, true);
    Spack sat_liq_p = Functions::polysvp1(temps, false);

    Spack mix_ice_r = Functions::qv_sat(temps, pres, true);
    Spack mix_liq_r = Functions::qv_sat(temps, pres, false);

    for(int s = 0; s < sat_ice_p.n; ++s){
      // Test vapor pressure
      if (abs(sat_ice_p[s] - correct_sat_ice_p) > C::Tol ) {errors++;}
      if (abs(sat_liq_p[s] - correct_sat_liq_p) > C::Tol) {errors++;}
      //Test mixing-ratios
      if (abs(mix_ice_r[s] -  correct_mix_ice_r) > C::Tol ) {errors++;}
      if (abs(mix_liq_r[s] -  correct_mix_liq_r) > C::Tol ) {errors++;}
    }
  }

  static void run()
  {
    int nerr = 0;
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("TestTableIce::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {

      errors = 0;

      // Test values @ the melting point of H20 @ 1e5 Pa
      saturation_tests(C::Tmelt, 1e5, 610.7960763188032, 610.7960763188032,
         0.003822318507864685,  0.003822318507864685, errors);

      //Test vaules @ 243.15K @ 1e5 Pa
      saturation_tests(243.15, 1e5, 37.98530141245404, 50.98455924912173,
         0.00023634717905493638,  0.0003172707211143376, errors);

      //Test values @ 303.15 @ 1e5 Pa
      saturation_tests(303.15, 1e5, 4242.757341329608, 4242.757341329608,
        0.0275579183092878, 0.0275579183092878, errors);

    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
};

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

TEST_CASE("p3_ice_tables", "[p3_functions]")
{
  UnitWrap::UnitTest<scream::DefaultDevice>::TestTableIce::test_read_lookup_tables_bfb();
  UnitWrap::UnitTest<scream::DefaultDevice>::TestTableIce::run_phys();
  UnitWrap::UnitTest<scream::DefaultDevice>::TestTableIce::run_bfb();
}

TEST_CASE("p3_functions", "[p3_functions]")
{
  UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Func::run();
}

} // namespace
