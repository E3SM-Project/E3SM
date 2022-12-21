#include "catch2/catch.hpp"

#include "physics/share/physics_test_data.hpp"
#include "share/scream_types.hpp"
#include "physics_unit_tests_common.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

namespace scream {
namespace physics {
namespace unit_test {

namespace {

struct FakeClass1 : public PhysicsTestData
{
  Int dim1, dim2, dim3, dim4;

  Real *one123, *two123, *three124;
  Real scalar;

  Int *int12, *int124;

  bool *bools12, *bools1;

  FakeClass1(Int dim1_, Int dim2_, Int dim3_, Int dim4_, Real scalar_) :
    PhysicsTestData({ {dim1_, dim2_, dim3_}, {dim1_, dim2_, dim4_}, {dim1_, dim2_},  {dim1_, dim2_, dim4_}, {dim1_, dim2_}, {dim1_}},
                    { {&one123, &two123},    {&three124},        }, {{&int12},       {&int124} },         {{&bools12},      {&bools1}}),
    dim1(dim1_),
    dim2(dim2_),
    dim3(dim3_),
    dim4(dim4_),
    scalar(scalar_)
  {}

  PTD_STD_DEF(FakeClass1, 5, dim1, dim2, dim3, dim4, scalar);
};

struct AllTrue : public PhysicsTestData
{
  Int dim1;
  bool *bools;

  AllTrue(Int dim1_) :
    PhysicsTestData({ {dim1_} }, {/*no reals*/}, {/*no ints*/},
                    { {&bools} }),
    dim1(dim1_)
  {}

  PTD_STD_DEF(AllTrue, 1, dim1);
};

} // empty namespace

template <typename D>
struct UnitWrap::UnitTest<D>::TestTestData
{
  static void test_fake_class1()
  {
    auto engine = setup_random_test();

    const std::vector<std::tuple<Int, Int, Int, Int>> dims = { std::make_tuple(7, 13, 29, 3), std::make_tuple(4, 8, 16, 2) };
    const std::vector<Real> scalars = { 42.1, 43.1 };
    FakeClass1 fakes_1[] = {
      FakeClass1(std::get<0>(dims[0]), std::get<1>(dims[0]), std::get<2>(dims[0]), std::get<3>(dims[0]), scalars[0]),
      FakeClass1(std::get<0>(dims[1]), std::get<1>(dims[1]), std::get<2>(dims[1]), std::get<3>(dims[1]), scalars[1]),
    };

    for (auto& d : fakes_1) {
      d.randomize(engine, { {d.two123, {-2.0, -1.0}}, {d.three124, {-3.0, -2.0}}, {d.int124, {42, 84}} });
    }

    static constexpr Int num_runs = sizeof(fakes_1) / sizeof(FakeClass1);

    // Copy construction
    FakeClass1 fakes_2[] = {
      FakeClass1(fakes_1[0]),
      FakeClass1(fakes_1[1]),
    };

    FakeClass1 fakes_3[] = {
      FakeClass1(1, 1, 1, 1, 1),
      FakeClass1(1, 1, 1, 1, 1),
    };

    // Assignment
    for (Int n = 0; n < num_runs; ++n) {
      fakes_3[n] = fakes_2[n];
    }

    for (Int n = 0; n < num_runs; ++n) {
      auto& d1 = fakes_1[n];
      auto& d2 = fakes_2[n];
      auto& d3 = fakes_3[n];

      // Check dimensions
      REQUIRE(d1.total(d1.one123) == std::get<0>(dims[n]) * std::get<1>(dims[n]) * std::get<2>(dims[n]));
      REQUIRE(d1.total(d1.one123) == d2.total(d2.one123));
      REQUIRE(d1.total(d1.one123) == d3.total(d3.one123));

      REQUIRE(d1.total(d1.two123) == std::get<0>(dims[n]) * std::get<1>(dims[n]) * std::get<2>(dims[n]));
      REQUIRE(d1.total(d1.two123) == d2.total(d2.two123));
      REQUIRE(d1.total(d1.two123) == d3.total(d3.two123));

      REQUIRE(d1.total(d1.three124) == std::get<0>(dims[n]) * std::get<1>(dims[n]) * std::get<3>(dims[n]));
      REQUIRE(d1.total(d1.three124) == d2.total(d2.three124));
      REQUIRE(d1.total(d1.three124) == d3.total(d3.three124));

      REQUIRE(d1.total(d1.int12) == std::get<0>(dims[n]) * std::get<1>(dims[n]));
      REQUIRE(d1.total(d1.int12) == d2.total(d2.int12));
      REQUIRE(d1.total(d1.int12) == d3.total(d3.int12));

      REQUIRE(d1.total(d1.int124) == std::get<0>(dims[n]) * std::get<1>(dims[n]) * std::get<3>(dims[n]));
      REQUIRE(d1.total(d1.int124) == d2.total(d2.int124));
      REQUIRE(d1.total(d1.int124) == d3.total(d3.int124));

      REQUIRE(d1.total(d1.bools12) == std::get<0>(dims[n]) * std::get<1>(dims[n]));
      REQUIRE(d1.total(d1.bools12) == d2.total(d2.bools12));
      REQUIRE(d1.total(d1.bools12) == d3.total(d3.bools12));
      REQUIRE(d1.total(d1.bools12) == d1.total(d1.int12));

      REQUIRE(d1.total(d1.bools1) == std::get<0>(dims[n]));
      REQUIRE(d1.total(d1.bools1) == d2.total(d2.bools1));
      REQUIRE(d1.total(d1.bools1) == d3.total(d3.bools1));

      // Check scalars

      REQUIRE(d1.dim1 == std::get<0>(dims[n]));
      REQUIRE(d1.dim1 == d2.dim1);
      REQUIRE(d1.dim1 == d3.dim1);

      REQUIRE(d1.dim2 == std::get<1>(dims[n]));
      REQUIRE(d1.dim2 == d2.dim2);
      REQUIRE(d1.dim2 == d3.dim2);

      REQUIRE(d1.dim3 == std::get<2>(dims[n]));
      REQUIRE(d1.dim3 == d2.dim3);
      REQUIRE(d1.dim3 == d3.dim3);

      REQUIRE(d1.dim4 == std::get<3>(dims[n]));
      REQUIRE(d1.dim4 == d2.dim4);
      REQUIRE(d1.dim4 == d3.dim4);

      REQUIRE(d1.scalar == scalars[n]);
      REQUIRE(d1.scalar == d2.scalar);
      REQUIRE(d1.scalar == d3.scalar);

      // Check randomization and correct copy construction, assignment
      for (Int i = 0; i < d1.total(d1.one123); ++i) {
        REQUIRE( (d1.one123[i] >= 0.0  && d1.one123[i] <= 1.0) );
        REQUIRE( (d1.two123[i] >= -2.0 && d1.two123[i] <= -1.0) );

        REQUIRE(d1.one123[i] == d2.one123[i]);
        REQUIRE(d1.one123[i] == d3.one123[i]);

        REQUIRE(d1.two123[i] == d2.two123[i]);
        REQUIRE(d1.two123[i] == d3.two123[i]);
      }
      for (Int i = 0; i < d1.total(d1.three124); ++i) {
        REQUIRE( (d1.three124[i] >= -3.0 && d1.three124[i] <= -2.0) );
        REQUIRE( (d1.int124[i]  >= 42   && d1.int124[i]  <= 84) );

        REQUIRE(d1.three124[i] == d2.three124[i]);
        REQUIRE(d1.three124[i] == d3.three124[i]);

        REQUIRE(d1.int124[i] == d2.int124[i]);
        REQUIRE(d1.int124[i] == d3.int124[i]);
      }
      for (Int i = 0; i < d1.total(d1.int12); ++i) {
        REQUIRE( (d1.int12[i] >= 0 && d1.int12[i] <= 1) );

        REQUIRE(d1.int12[i] == d2.int12[i]);
        REQUIRE(d1.int12[i] == d3.int12[i]);

        REQUIRE( (d1.bools12[i] == false || d1.bools12[i] == true) );

        REQUIRE(d1.bools12[i] == d2.bools12[i]);
        REQUIRE(d1.bools12[i] == d3.bools12[i]);
      }
      for (Int i = 0; i < d1.dim1; ++i) {
        REQUIRE( (d1.bools1[i] == false || d1.bools1[i] == true) );

        REQUIRE(d1.bools1[i] == d2.bools1[i]);
        REQUIRE(d1.bools1[i] == d3.bools1[i]);
      }

      // Check transpose
      d1.transpose<ekat::TransposeDirection::c2f>();
      for (Int i = 0; i < d1.dim(d1.one123, 0); ++i) {
        for (Int j = 0; j < d1.dim(d1.one123, 1); ++j) {
          const Int cidx = d1.dim2*i + j;
          const Int fidx = d1.dim1*j + i;

          REQUIRE(d1.int12[fidx] == d2.int12[cidx]+1);

          for (Int k = 0; k < d1.dim(d1.one123, 2); ++k) {
            const Int c3didx = (d1.dim2*d1.dim3)*i + j*d1.dim3 + k;
            const Int f3didx = (d1.dim1*d1.dim2)*k + j*d1.dim1 + i;
            REQUIRE(d1.one123[f3didx] == d2.one123[c3didx]);
            REQUIRE(d1.two123[f3didx] == d2.two123[c3didx]);
          }
          for (Int k = 0; k < d1.dim(d1.three124, 2); ++k) {
            const Int c3didx = (d1.dim2*d1.dim4)*i + j*d1.dim4 + k;
            const Int f3didx = (d1.dim1*d1.dim2)*k + j*d1.dim1 + i;
            REQUIRE(d1.three124[f3didx] == d2.three124[c3didx]);
          }
        }
      }

      d1.transpose<ekat::TransposeDirection::f2c>();
      for (Int i = 0; i < d1.total(d1.one123); ++i) {
        REQUIRE(d1.one123[i] == d2.one123[i]);
        REQUIRE(d1.two123[i] == d2.two123[i]);
      }
      for (Int i = 0; i < d1.total(d1.three124); ++i) {
        REQUIRE(d1.three124[i] == d2.three124[i]);
      }
      for (Int i = 0; i < d1.total(d1.int12); ++i) {
        REQUIRE(d1.int12[i] == d2.int12[i]);
      }

      // Check that different data obj are using different memory spaces
      d1.one123[0] = 123.0;
      d1.int12[0]  = 123;

      REQUIRE(d1.one123[0] != d2.one123[0]);
      REQUIRE(d1.one123[0] != d3.one123[0]);

      REQUIRE(d1.int12[0] != d2.int12[0]);
      REQUIRE(d1.int12[0] != d3.int12[0]);
    }
  }

  static void test_all_true()
  {
    auto engine = setup_random_test();

    static constexpr int size = 42;
    AllTrue at(size);

    at.randomize(engine, { {at.bools, {1.0, 1.0}} });

    for (int i = 0; i < size; ++i) {
      REQUIRE(at.bools[i]);
    }
  }

  static void run()
  {
    test_fake_class1();
    test_all_true();
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
