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
  Real *one12, *two12, *three13, *four13, *five3d;

  Int* ints;

  FakeClass1(Int dim1, Int dim2, Int dim3) :
    PhysicsTestData(dim1, dim2, dim3,
                    {&one12, &two12}, {&three13, &four13}, {}, {&ints}, {&five3d}) {}

  PTD_STD_DEF(FakeClass1, 3, 0);
};

struct FakeClass2 : public PhysicsTestData
{
  Real *one12, *two1;

  FakeClass2(Int dim1, Int dim2) :
    PhysicsTestData(dim1, dim2, {&one12}, {&two1}) {}

  PTD_STD_DEF(FakeClass2, 2, 0);
};

struct FakeClass3 : public PhysicsTestData
{
  Real *one1, *two1;
  Real scalar;

  FakeClass3(Int dim1, Real scalar_) :
    PhysicsTestData(dim1, {&one1, &two1}), scalar(scalar_) {}

  PTD_DIM_RENAME(1, foo);
  PTD_STD_DEF(FakeClass3, 1, 1, scalar);
};

struct FakeClass4 : public PhysicsTestDataGeneric
{
  Int dim1, dim2, dim3, dim4;

  Real *one123, *two123, *three124;
  Real scalar;

  Int *int12, *int124;

  FakeClass4(Int dim1_, Int dim2_, Int dim3_, Int dim4_, Real scalar_) :
    PhysicsTestDataGeneric({ {dim1_, dim2_, dim3_}, {dim1_, dim2_, dim4_}, {dim1_, dim2_},  {dim1_, dim2_, dim4_} },
                           { {&one123, &two123},    {&three124},        }, {{&int12},       {&int124} }),
    dim1(dim1_),
    dim2(dim2_),
    dim3(dim3_),
    dim4(dim4_),
    scalar(scalar_)
  {}

  PTDG_STD_DEF(FakeClass4, 5, dim1, dim2, dim3, dim4, scalar);
};

} // empty namespace

template <typename D>
struct UnitWrap::UnitTest<D>::TestTestData
{
  static void test_fake_class1()
  {
    const std::vector<std::tuple<Int, Int, Int>> dims = { std::make_tuple(7, 13, 29), std::make_tuple(4, 8, 16) };
    FakeClass1 fakes1_1[] = {
      FakeClass1(std::get<0>(dims[0]), std::get<1>(dims[0]), std::get<2>(dims[0])),
      FakeClass1(std::get<0>(dims[1]), std::get<1>(dims[1]), std::get<2>(dims[1])),
    };

    for (auto& d : fakes1_1) {
      d.randomize({ {d.two12, {-2.0, -1.0}}, {d.three13, {-3.0, -2.0}}, {d.ints, {42, 84}}, {d.five3d, {6.0, 9.0}} });
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
      REQUIRE(d1.total1x2() == std::get<0>(dims[n]) * std::get<1>(dims[n]));
      REQUIRE(d1.total1x2() == d2.total1x2());
      REQUIRE(d1.total1x2() == d3.total1x2());

      REQUIRE(d1.total1x3() == std::get<0>(dims[n]) * std::get<2>(dims[n]));
      REQUIRE(d1.total1x3() == d2.total1x3());
      REQUIRE(d1.total1x3() == d3.total1x3());

      REQUIRE(d1.total3d() == std::get<0>(dims[n]) * std::get<1>(dims[n]) * std::get<2>(dims[n]));
      REQUIRE(d1.total3d() == d2.total3d());
      REQUIRE(d1.total3d() == d3.total3d());

      REQUIRE(d1.dim1 == std::get<0>(dims[n]));
      REQUIRE(d1.dim1 == d2.dim1);
      REQUIRE(d1.dim1 == d3.dim1);

      // Check randomization and correct copy construction, assignment
      for (Int i = 0; i < d1.total1x2(); ++i) {
        REQUIRE( (d1.one12[i] >= 0.0  && d1.one12[i] <= 1.0) );
        REQUIRE( (d1.two12[i] >= -2.0 && d1.two12[i] <= -1.0) );

        REQUIRE(d1.one12[i] == d2.one12[i]);
        REQUIRE(d1.one12[i] == d3.one12[i]);

        REQUIRE(d1.two12[i] == d2.two12[i]);
        REQUIRE(d1.two12[i] == d3.two12[i]);
      }
      for (Int i = 0; i < d1.total1x3(); ++i) {
        REQUIRE( (d1.three13[i] >= -3.0 && d1.three13[i] <= -2.0) );
        REQUIRE( (d1.four13[i] >= 0.0   && d1.four13[i] <= 1.0) );

        REQUIRE(d1.three13[i] == d2.three13[i]);
        REQUIRE(d1.three13[i] == d3.three13[i]);

        REQUIRE(d1.four13[i] == d2.four13[i]);
        REQUIRE(d1.four13[i] == d3.four13[i]);
      }
      for (Int i = 0; i < d1.total3d(); ++i) {
        REQUIRE( (d1.five3d[i] >= 6.0 && d1.five3d[i] <= 9.0) );

        REQUIRE(d1.five3d[i] == d2.five3d[i]);
        REQUIRE(d1.five3d[i] == d3.five3d[i]);
      }
      for (Int i = 0; i < d1.dim1; ++i) {
        REQUIRE( (d1.ints[i] >= 42 && d1.ints[i] <= 84) );

        REQUIRE(d1.ints[i] == d2.ints[i]);
        REQUIRE(d1.ints[i] == d3.ints[i]);
      }

      // Check transpose
      d1.transpose<ekat::TransposeDirection::c2f>();
      for (Int i = 0; i < d1.dim1; ++i) {
        for (Int j = 0; j < d1.dim2; ++j) {
          const Int cidx = d1.dim2*i + j;
          const Int fidx = d1.dim1*j + i;

          REQUIRE(d1.one12[fidx] == d2.one12[cidx]);
          REQUIRE(d1.two12[fidx] == d2.two12[cidx]);

          for (Int k = 0; k < d1.dim3; ++k) {
            const Int c3didx = (d1.dim2*d1.dim3)*i + j*d1.dim3 + k;
            const Int f3didx = (d1.dim1*d1.dim2)*k + j*d1.dim1 + i;
            REQUIRE(d1.five3d[f3didx] == d2.five3d[c3didx]);
          }
        }
        for (Int j = 0; j < d1.dim3; ++j) {
          const Int cidx = d1.dim3*i + j;
          const Int fidx = d1.dim1*j + i;

          REQUIRE(d1.three13[fidx] == d2.three13[cidx]);
          REQUIRE(d1.four13[fidx] == d2.four13[cidx]);
        }
      }

      d1.transpose<ekat::TransposeDirection::f2c>();
      for (Int i = 0; i < d1.total1x2(); ++i) {
        REQUIRE(d1.one12[i] == d2.one12[i]);
        REQUIRE(d1.two12[i] == d2.two12[i]);
      }
      for (Int i = 0; i < d1.total1x3(); ++i) {
        REQUIRE(d1.three13[i] == d2.three13[i]);
        REQUIRE(d1.four13[i]  == d2.four13[i]);
      }
      for (Int i = 0; i < d1.total3d(); ++i) {
        REQUIRE(d1.five3d[i] == d2.five3d[i]);
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

  static void test_fake_class2()
  {
    const std::vector<std::tuple<Int, Int>> dims = { std::make_tuple(7, 13), std::make_tuple(4, 8) };
    FakeClass2 fakes1_1[] = {
      FakeClass2(std::get<0>(dims[0]), std::get<1>(dims[0])),
      FakeClass2(std::get<0>(dims[1]), std::get<1>(dims[1])),
    };

    for (auto& d : fakes1_1) {
      d.randomize({ {d.two1, {-2.0, -1.0}} });
    }

    static constexpr Int num_runs = sizeof(fakes1_1) / sizeof(FakeClass2);

    // Copy construction
    FakeClass2 fakes1_2[] = {
      FakeClass2(fakes1_1[0]),
      FakeClass2(fakes1_1[1]),
    };

    FakeClass2 fakes1_3[] = {
      FakeClass2(1, 1),
      FakeClass2(1, 1),
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
      REQUIRE(d1.total1x2() == std::get<0>(dims[n]) * std::get<1>(dims[n]));
      REQUIRE(d1.total1x2() == d2.total1x2());
      REQUIRE(d1.total1x2() == d3.total1x2());

      REQUIRE(d1.dim1 == std::get<0>(dims[n]));
      REQUIRE(d1.dim1 == d2.dim1);
      REQUIRE(d1.dim1 == d3.dim1);

      // Check randomization and correct copy construction, assignment
      for (Int i = 0; i < d1.total1x2(); ++i) {
        REQUIRE( (d1.one12[i] >= 0.0 && d1.one12[i] <= 1.0) );

        REQUIRE(d1.one12[i] == d2.one12[i]);
        REQUIRE(d1.one12[i] == d3.one12[i]);
      }
      for (Int i = 0; i < d1.dim1; ++i) {
        REQUIRE( (d1.two1[i] >= -2.0 && d1.two1[i] <= -1.0) );

        REQUIRE(d1.two1[i] == d2.two1[i]);
        REQUIRE(d1.two1[i] == d3.two1[i]);
      }

      // Check transpose
      d1.transpose<ekat::TransposeDirection::c2f>();
      for (Int i = 0; i < d1.dim1; ++i) {
        for (Int j = 0; j < d1.dim2; ++j) {
          const Int cidx = d1.dim2*i + j;
          const Int fidx = d1.dim1*j + i;

          REQUIRE(d1.one12[fidx] == d2.one12[cidx]);
        }
      }

      d1.transpose<ekat::TransposeDirection::f2c>();
      for (Int i = 0; i < d1.total1x2(); ++i) {
        REQUIRE(d1.one12[i] == d2.one12[i]);
      }

      // Check that different data obj are using different memory spaces
      d1.one12[0] = 123.0;
      d1.two1[0]  = 123;

      REQUIRE(d1.one12[0] != d2.one12[0]);
      REQUIRE(d1.one12[0] != d3.one12[0]);

      REQUIRE(d1.two1[0] != d2.two1[0]);
      REQUIRE(d1.two1[0] != d3.two1[0]);
    }
  }

  static void test_fake_class3()
  {
    const std::vector<Int> dims = { 7, 8 };
    const std::vector<Real> scalars = { 42.1, 43.1 };
    FakeClass3 fakes1_1[] = {
      FakeClass3(dims[0], scalars[0]),
      FakeClass3(dims[1], scalars[1]),
    };

    for (auto& d : fakes1_1) {
      d.randomize({ {d.two1, {-2.0, -1.0}} });
    }

    static constexpr Int num_runs = sizeof(fakes1_1) / sizeof(FakeClass3);

    // Copy construction
    FakeClass3 fakes1_2[] = {
      FakeClass3(fakes1_1[0]),
      FakeClass3(fakes1_1[1]),
    };

    FakeClass3 fakes1_3[] = {
      FakeClass3(1, 0),
      FakeClass3(1, 0),
    };

    // Assignment
    for (Int n = 0; n < num_runs; ++n) {
      fakes1_3[n] = fakes1_2[n];
    }

    for (Int n = 0; n < num_runs; ++n) {
      auto& d1 = fakes1_1[n];
      auto& d2 = fakes1_2[n];
      auto& d3 = fakes1_3[n];

      // Check dim rename
      REQUIRE(d1.foo() == d1.dim1);

      // Check dimensions
      REQUIRE(d1.dim1 == dims[n]);
      REQUIRE(d1.dim1 == d2.dim1);
      REQUIRE(d1.dim1 == d3.dim1);

      // Check scalar
      REQUIRE(d1.scalar == scalars[n]);
      REQUIRE(d1.scalar == d2.scalar);
      REQUIRE(d1.scalar == d3.scalar);

      // Check randomization and correct copy construction, assignment
      for (Int i = 0; i < d1.dim1; ++i) {
        REQUIRE( (d1.one1[i] >= 0.0  && d2.one1[i] <= 1.0) );
        REQUIRE( (d1.two1[i] >= -2.0 && d1.two1[i] <= -1.0) );
      }

      // Check that different data obj are using different memory spaces
      d1.one1[0] = 123.0;

      REQUIRE(d1.one1[0] != d2.one1[0]);
      REQUIRE(d1.one1[0] != d3.one1[0]);
    }
  }

  static void test_fake_class4()
  {
    const std::vector<std::tuple<Int, Int, Int, Int>> dims = { std::make_tuple(7, 13, 29, 3), std::make_tuple(4, 8, 16, 2) };
    const std::vector<Real> scalars = { 42.1, 43.1 };
    FakeClass4 fakes_1[] = {
      FakeClass4(std::get<0>(dims[0]), std::get<1>(dims[0]), std::get<2>(dims[0]), std::get<3>(dims[0]), scalars[0]),
      FakeClass4(std::get<0>(dims[1]), std::get<1>(dims[1]), std::get<2>(dims[1]), std::get<3>(dims[1]), scalars[1]),
    };

    for (auto& d : fakes_1) {
      d.randomize({ {d.two123, {-2.0, -1.0}}, {d.three124, {-3.0, -2.0}}, {d.int124, {42, 84}} });
    }

    static constexpr Int num_runs = sizeof(fakes_1) / sizeof(FakeClass1);

    // Copy construction
    FakeClass4 fakes_2[] = {
      FakeClass4(fakes_1[0]),
      FakeClass4(fakes_1[1]),
    };

    FakeClass4 fakes_3[] = {
      FakeClass4(1, 1, 1, 1, 1),
      FakeClass4(1, 1, 1, 1, 1),
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

  static void run()
  {
    test_fake_class1();
    test_fake_class2();
    test_fake_class3();
    test_fake_class4();
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
