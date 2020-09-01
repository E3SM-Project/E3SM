#include "catch2/catch.hpp"

#include "physics/share/physics_test_data.hpp"
#include "share/scream_types.hpp"
#include "physics_unit_tests_common.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace scream {
namespace physics {
namespace unit_test {

namespace {

struct FakeClass1 : public PhysicsTestData
{
  Real *one12, *two12, *three13, *four13;

  Int* ints;

  FakeClass1(Int dim1, Int dim2, Int dim3) :
    PhysicsTestData(dim1, dim2, dim3,
                    {&one12, &two12}, {&three13, &four13}, {}, {&ints}) {}

  FakeClass1(const FakeClass1& rhs) :
    PhysicsTestData(rhs, {&one12, &two12, &three13, &four13}, {&ints}) {}

  FakeClass1& operator=(const FakeClass1& rhs) { PhysicsTestData::operator=(rhs); return *this; }
};

} // empty namespace

template <typename D>
struct UnitWrap::UnitTest<D>::TestTestData
{
  static void run()
  {
    const std::vector<std::tuple<Int, Int, Int>> dims = { std::make_tuple(7, 13, 29), std::make_tuple(4, 8, 16) };
    FakeClass1 fakes1_1[] = {
      FakeClass1(std::get<0>(dims[0]), std::get<1>(dims[0]), std::get<2>(dims[0])),
      FakeClass1(std::get<0>(dims[1]), std::get<1>(dims[1]), std::get<2>(dims[1])),
    };

    for (auto& d : fakes1_1) {
      d.randomize({ {d.two12, {-2.0, -1.0}}, {d.three13, {-3.0, -2.0}}, {d.ints, {42, 84}} });
    }

    static constexpr Int num_runs = sizeof(fakes1_1) / sizeof(FakeClass1);

    // Copy construction
    FakeClass1 fakes1_2[] = {
      FakeClass1(fakes1_1[0]),
      FakeClass1(fakes1_1[1]),
    };

    FakeClass1 fakes1_3[] = {
      FakeClass1(1, 1, 1),
      FakeClass1(1, 1, 1),
    };

    // Assignment
    for (Int n = 0; n < num_runs; ++n) {
      fakes1_3[n] = fakes1_2[n];
    }

    for (Int n = 0; n < num_runs; ++n) {
      auto& d1 = fakes1_1[n];
      auto& d2 = fakes1_2[n];
      auto& d3 = fakes1_3[n];

      // Check dimensions
      REQUIRE(d1.total() == std::get<0>(dims[n]) * std::get<1>(dims[n]));
      REQUIRE(d1.total() == d2.total());
      REQUIRE(d1.total() == d3.total());

      REQUIRE(d1.totali() == std::get<0>(dims[n]) * std::get<2>(dims[n]));
      REQUIRE(d1.totali() == d2.totali());
      REQUIRE(d1.totali() == d3.totali());

      REQUIRE(d1.shcol == std::get<0>(dims[n]));
      REQUIRE(d1.shcol == d2.shcol);
      REQUIRE(d1.shcol == d3.shcol);

      // Check randomization and correct copy construction, assignment
      for (Int i = 0; i < d1.total(); ++i) {
        REQUIRE( (d1.one12[i] > 0.0  && d1.one12[i] < 1.0) );
        REQUIRE( (d1.two12[i] > -2.0 && d1.two12[i] < -1.0) );

        REQUIRE(d1.one12[i] == d2.one12[i]);
        REQUIRE(d1.one12[i] == d3.one12[i]);

        REQUIRE(d1.two12[i] == d2.two12[i]);
        REQUIRE(d1.two12[i] == d3.two12[i]);
      }
      for (Int i = 0; i < d1.totali(); ++i) {
        REQUIRE( (d1.three13[i] > -3.0 && d1.three13[i] < -2.0) );
        REQUIRE( (d1.four13[i] > 0.0   && d1.four13[i] < 1.0) );

        REQUIRE(d1.three13[i] == d2.three13[i]);
        REQUIRE(d1.three13[i] == d3.three13[i]);

        REQUIRE(d1.four13[i] == d2.four13[i]);
        REQUIRE(d1.four13[i] == d3.four13[i]);
      }
      for (Int i = 0; i < d1.shcol; ++i) {
        REQUIRE( (d1.ints[i] >= 42 && d1.ints[i] <= 84) );

        REQUIRE(d1.ints[i] == d2.ints[i]);
        REQUIRE(d1.ints[i] == d3.ints[i]);
      }

      // Check transpose
      d1.transpose<ekat::util::TransposeDirection::c2f>();
      for (Int i = 0; i < d1.shcol; ++i) {
        for (Int j = 0; j < d1.nlev; ++j) {
          const Int cidx = d1.nlev*i + j;
          const Int fidx = d1.shcol*j + i;

          REQUIRE(d1.one12[fidx] == d2.one12[cidx]);
          REQUIRE(d1.two12[fidx] == d2.two12[cidx]);
        }
        for (Int j = 0; j < d1.nlevi; ++j) {
          const Int cidx = d1.nlevi*i + j;
          const Int fidx = d1.shcol*j + i;

          REQUIRE(d1.three13[fidx] == d2.three13[cidx]);
          REQUIRE(d1.four13[fidx] == d2.four13[cidx]);
        }
      }

      d1.transpose<ekat::util::TransposeDirection::f2c>();
      for (Int i = 0; i < d1.total(); ++i) {
        REQUIRE(d1.one12[i] == d2.one12[i]);
        REQUIRE(d1.two12[i] == d2.two12[i]);
      }
      for (Int i = 0; i < d1.totali(); ++i) {
        REQUIRE(d1.three13[i] == d2.three13[i]);
        REQUIRE(d1.four13[i]  == d2.four13[i]);
      }

      // Check that different data obj are using different memory spaces
      d1.one12[0] = 123.0;
      d1.ints[0]  = 123;

      REQUIRE(d1.one12[0] != d2.one12[0]);
      REQUIRE(d1.one12[0] != d3.one12[0]);

      REQUIRE(d1.ints[0] != d2.ints[0]);
      REQUIRE(d1.ints[0] != d3.ints[0]);
    }
  }

};

} // namespace unit_test
} // namespace physics
} // namespace scream

namespace {

TEST_CASE("physics_test_data", "[physics_test_data]")
{
  scream::physics::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestTestData::run();
}

} // namespace
