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

/*
 * Unit-tests for p3 ice table functions.
 */

template <typename D>
struct UnitWrap::UnitTest<D>::TestTableIce {

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
    using KTH = KokkosTypes<HostDevice>;

    // Read in ice tables
    view_itab_table itab;
    view_itabcol_table itabcol;
    Functions::init_kokkos_ice_lookup_tables(itab, itabcol);

    constexpr Scalar qsmall = C::QSMALL;
    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    // Load some lookup inputs, need at least one per pack value
    LookupIceData lid[max_pack_size] = {
      // qitot,   nitot,     qirim,     rhop
      {0.971E-07, 0.657E+06, 0.971E-07, 0.900E+03},
      {0.510E-02, 0.454E+06, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+07, 0.811E-10, 0.500E+02},

      {0.971E-07, 0.657E+06, 0.971E-07, 0.900E+03},
      {0.510E-02, 0.454E+01, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+06, 0.811E-10, 0.500E+02},

      {0.971E-07, 0.657E+06, 0.271E-06, 0.900E+03},
      {0.510E-02, 0.454E+06, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+06, 0.811E-10, 0.500E+02},

      {0.971E-07, 0.657E+06, 0.971E-07, 0.200E+04},
      {0.510E-02, 0.454E+06, 0.714E-05, 0.500E+02},
      {0.500E-07, 0.545E+06, 0.000E+00, 0.000E+00},
      {0.136E-08, 0.487E+06, 0.811E-10, 0.500E+02}
    };

    LookupIceDataB lidb[max_pack_size] = {
      // qr,      nr
      {0.263E-05, 0.100E+07},
      {0.100E-01, 0.100E+07},
      {0.000E+00, 0.100E-15},
      {0.263E-05, 0.100E+07},

      {0.263E-05, 0.100E+07},
      {0.100E-01, 0.100E+07},
      {0.000E+00, 0.0      },
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
    for (Int i = 0; i < Spack::n; ++i) {
      find_lookuptable_indices_1a(lid[i]);
      find_lookuptable_indices_1b(lidb[i]);
      access_lookup_table(altd[i]);
      access_lookup_table_coll(altcd[i]);
    }

    // Sync to device
    KTH::view_1d<LookupIceData> lid_host("lid_host", Spack::n);
    KTH::view_1d<LookupIceDataB> lidb_host("lidb_host", Spack::n);
    view_1d<LookupIceData> lid_device("lid_device", Spack::n);
    view_1d<LookupIceDataB> lidb_device("lidb_device", Spack::n);
    std::copy(&lid[0], &lid[0] + Spack::n, lid_host.data());
    std::copy(&lidb[0], &lidb[0] + Spack::n, lidb_host.data());
    Kokkos::deep_copy(lid_device, lid_host);
    Kokkos::deep_copy(lidb_device, lidb_host);

    // Run the lookup from a kernel and copy results back to host
    view_1d<IntSmallPack> int_results("int results", 5);
    view_1d<Spack> real_results("real results", 7);
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      // Init packs
      TableIce ti;
      TableRain tr;
      Spack qitot, nitot, qirim, rhop, qr, nr;
      for(Int s = 0; s < Spack::n; ++s) {
        qitot[s] = lid_device(s).qitot;
        nitot[s] = lid_device(s).nitot;
        qirim[s] = lid_device(s).qirim;
        rhop[s]  = lid_device(s).rhop;

        qr[s]    = lidb_device(s).qr;
        nr[s]    = lidb_device(s).nr;
      }

      Smask qiti_gt_small(qitot > qsmall);
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

    // Sync results back to host
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
#if 0
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
#endif
  }
};

}
}
}

namespace {

TEST_CASE("p3_ice_tables", "[p3_functions]")
{
  using TTI = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestTableIce;

  TTI::test_read_lookup_tables_bfb();
  TTI::run_phys();
  TTI::run_bfb();
}

}
