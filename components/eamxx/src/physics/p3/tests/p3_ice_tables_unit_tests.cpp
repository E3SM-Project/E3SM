#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

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
struct UnitWrap::UnitTest<D>::TestTableIce : public UnitWrap::UnitTest<D>::Base {

  void test_read_lookup_tables_bfb()
  {
    // Read in ice tables
    view_ice_table ice_table_vals;
    view_collect_table collect_table_vals;
    Functions::init_kokkos_ice_lookup_tables(ice_table_vals, collect_table_vals);

    // Get data from fortran
    P3InitAP3Data d;
    p3_init_a(d);

    // Copy device data to host
    const auto ice_table_vals_host = Kokkos::create_mirror_view(ice_table_vals);
    const auto collect_table_vals_host = Kokkos::create_mirror_view(collect_table_vals);
    Kokkos::deep_copy(ice_table_vals_host, ice_table_vals);
    Kokkos::deep_copy(collect_table_vals_host, collect_table_vals);

    // Compare (on host)
    for (size_t i = 0; i < ice_table_vals_host.extent(0); ++i) {
      for (size_t j = 0; j < ice_table_vals_host.extent(1); ++j) {
        for (size_t k = 0; k < ice_table_vals_host.extent(2); ++k) {

          for (size_t l = 0; l < ice_table_vals_host.extent(3); ++l) {
            REQUIRE(ice_table_vals_host(i, j, k, l) == d.ice_table_vals(i, j, k, l));
          }

          for (size_t l = 0; l < collect_table_vals_host.extent(3); ++l) {
            for (size_t m = 0; m < collect_table_vals_host.extent(4); ++m) {
              REQUIRE(collect_table_vals_host(i, j, k, l, m) == d.collect_table_vals(i, j, k, l, m));
            }
          }

        }
      }
    }
  }

  template <typename View>
  void init_table_linear_dimension(View& table, int linear_dimension)
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

  void run_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    // Read in ice tables
    view_ice_table ice_table_vals;
    view_collect_table collect_table_vals;
    Functions::init_kokkos_ice_lookup_tables(ice_table_vals, collect_table_vals);

    constexpr Scalar qsmall = C::QSMALL;

    // Load some lookup inputs, need at least one per pack value
    LookupIceData lid[max_pack_size] = {
      // qi,   ni,     qm,     rhop
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

    // Sync to device
    KTH::view_1d<LookupIceData> lid_host("lid_host", max_pack_size);
    KTH::view_1d<LookupIceDataB> lidb_host("lidb_host", max_pack_size);
    view_1d<LookupIceData> lid_device("lid_device", max_pack_size);
    view_1d<LookupIceDataB> lidb_device("lidb_device", max_pack_size);
    std::copy(&lid[0], &lid[0] + max_pack_size, lid_host.data());
    std::copy(&lidb[0], &lidb[0] + max_pack_size, lidb_host.data());
    Kokkos::deep_copy(lid_device, lid_host);
    Kokkos::deep_copy(lidb_device, lidb_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        lid[i].read(Base::m_fid);
        lidb[i].read(Base::m_fid);
        altd[i].read(Base::m_fid);
        altcd[i].read(Base::m_fid);
      }
    }

    // Run the lookup from a kernel and copy results back to host
    view_2d<Int>  int_results("int results", 5, max_pack_size);
    view_2d<Real> real_results("real results", 7, max_pack_size);
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init packs
      TableIce ti;
      TableRain tr;
      Spack qi, ni, qm, rhop, qr, nr;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qi[s]    = lid_device(vs).qi;
        ni[s]    = lid_device(vs).ni;
        qm[s]    = lid_device(vs).qm;
        rhop[s]  = lid_device(vs).rhop;
        qr[s]    = lidb_device(vs).qr;
        nr[s]    = lidb_device(vs).nr;
      }

      Smask qiti_gt_small(qi > qsmall);
      Functions::lookup_ice(qi, ni, qm, rhop, ti, qiti_gt_small);
      Functions::lookup_rain(qr, nr, tr, qiti_gt_small);
      Spack ice_result = Functions::apply_table_ice(access_table_index-1, ice_table_vals, ti, qiti_gt_small);
      Spack rain_result = Functions::apply_table_coll(access_table_index-1, collect_table_vals, ti, tr, qiti_gt_small);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        int_results(0, vs) = ti.dumi[s];
        int_results(1, vs) = ti.dumjj[s];
        int_results(2, vs) = ti.dumii[s];
        int_results(3, vs) = ti.dumzz[s];

        int_results(4, vs) = tr.dumj[s];

        real_results(0, vs) = ti.dum1[s];
        real_results(1, vs) = ti.dum4[s];
        real_results(2, vs) = ti.dum5[s];
        real_results(3, vs) = ti.dum6[s];

        real_results(4, vs) = tr.dum3[s];

        real_results(5, vs) = ice_result[s];

        real_results(6, vs) = rain_result[s];
      }
    });

    // Sync results back to host
    auto int_results_mirror  = Kokkos::create_mirror_view(int_results);
    auto real_results_mirror = Kokkos::create_mirror_view(real_results);
    Kokkos::deep_copy(int_results_mirror, int_results);
    Kokkos::deep_copy(real_results_mirror, real_results);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for(int s = 0; s < max_pack_size; ++s) {
        REQUIRE(int_results_mirror(0, s) == lid[s].dumi);
        REQUIRE(int_results_mirror(1, s) == lid[s].dumjj);
        REQUIRE(int_results_mirror(2, s) == lid[s].dumii);
        REQUIRE(int_results_mirror(3, s) == lid[s].dumzz);

        REQUIRE(int_results_mirror(4, s) == lidb[s].dumj);

        REQUIRE(real_results_mirror(0, s) == lid[s].dum1);
        REQUIRE(real_results_mirror(1, s) == lid[s].dum4);
        REQUIRE(real_results_mirror(2, s) == lid[s].dum5);
        REQUIRE(real_results_mirror(3, s) == lid[s].dum6);

        REQUIRE(real_results_mirror(4, s) == lidb[s].dum3);

        REQUIRE(real_results_mirror(5, s) == altd[s].proc);

        REQUIRE(real_results_mirror(6, s) == altcd[s].proc);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        lid[s].dumi = int_results_mirror(0, s);
        lid[s].dumjj = int_results_mirror(1, s);
        lid[s].dumii = int_results_mirror(2, s);
        lid[s].dumzz = int_results_mirror(3, s);

        lidb[s].dumj = int_results_mirror(4, s);

        lid[s].dum1 = real_results_mirror(0, s);
        lid[s].dum4 = real_results_mirror(1, s);
        lid[s].dum5 = real_results_mirror(2, s);
        lid[s].dum6 = real_results_mirror(3, s);

        lidb[s].dum3 = real_results_mirror(4, s);

        altd[s].proc = real_results_mirror(5, s);

        altcd[s].proc = real_results_mirror(6, s);

        lid[s].write(Base::m_fid);
        lidb[s].write(Base::m_fid);
        altd[s].write(Base::m_fid);
        altcd[s].write(Base::m_fid);
      }
    }
  }

  void run_phys()
  {
#if 0
    view_ice_table ice_table_vals;
    init_table_linear_dimension(ice_table_vals, 0);

    int nerr = 0;
    TeamPolicy policy(ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ice_table_vals.extent(0), ice_table_vals.extent(1)));
    Kokkos::parallel_reduce("TestTableIce::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {
      //int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, ice_table_vals.extent(1)), [&] (const int& j) {

        for (size_t k = 0; k < ice_table_vals.extent(2); ++k) {
          for (size_t l = 0; l < ice_table_vals.extent(3); ++l) {
            // Init packs to same value, TODO: how to pick use values?
            Spack qi(0.1), ni(0.2), qm(0.3), rhop(0.4), qr(0.5), nr(0.6);

            TableIce ti;
            TableRain tr;
            Functions::lookup_ice(qi, ni, qm, rhop, ti);
            Functions::lookup_rain(qr, nr, tr);

            /*Spack proc1 = */ Functions::apply_table_ice(1, ice_table_vals, ti);
            //Spack proc2 = Functions::apply_table_coll(1, collect_table_vals, ti, tr);

            // TODO: how to test?
          }
        }
      });
      errors = 0;

    }, nerr);

    view_collect_table collect_table_vals;

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
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestTableIce;

  T t;
  t.test_read_lookup_tables_bfb();
  t.run_phys();
  t.run_bfb();
}

}
