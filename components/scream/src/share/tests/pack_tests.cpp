#include "catch2/catch.hpp"

#include "share/scream_config.hpp"
#include "share/scream_types.hpp"
#include "share/scream_pack.hpp"

namespace {

template <int PACKN>
struct TestMask {
  typedef scream::pack::Mask<PACKN> Mask;
  typedef scream::pack::Pack<int, PACKN> Pack;

  static int sum_true (const Mask& m) {
    int sum1 = 0, sum2 = 0, sum3 = 0;
    scream_masked_loop(m) ++sum1;
    scream_masked_loop_no_force_vec(m) ++sum2;
    scream_masked_loop_no_vec(m) ++sum3;
    REQUIRE(sum1 == sum2);
    REQUIRE(sum2 == sum3);
    return sum1;
  }

  static void run () {
    {
      Mask m(false);
      REQUIRE( ! m.any());
    }
    {
      Mask m(true);
      REQUIRE(m.any());
      for (int i = 0; i < Mask::n; ++i) REQUIRE(m[i]);
      REQUIRE(sum_true(m) == Mask::n);
    }
    for (int i = 0; i < Mask::n; ++i) {
      Mask m(false);
      m.set(i, true);
      REQUIRE(sum_true(m) == 1);
    }
    {
      Pack a, b;
      for (int i = 0; i < Mask::n; ++i)
        a[i] = i + (i % 2);
      for (int i = 0; i < Mask::n; ++i)
        b[i] = i;
      const auto m1 = a > b;
      scream_masked_loop_no_vec(m1) REQUIRE(m1[s] == (s % 2 == 1));
      const auto m2 = b < a;
      scream_masked_loop_no_vec(m2) REQUIRE(m2[s] == (s % 2 == 1));
      REQUIRE(sum_true(m1 && m2) == Mask::n / 2);
      REQUIRE(sum_true(m1 || m2) == Mask::n / 2);
      REQUIRE(sum_true(m1 && ! m2) == 0);
      REQUIRE(sum_true(m1 || ! m2) == Mask::n);
    }
  }
};

TEST_CASE("Mask", "scream::pack") {
  TestMask<1>::run();
  TestMask<2>::run();
  TestMask<3>::run();
  TestMask<4>::run();
  TestMask<8>::run();
  TestMask<16>::run();
  TestMask<32>::run();
}

template <int PACKN>
struct TestPack {
  typedef scream::pack::Mask<PACKN> Mask;
  typedef scream::pack::Pack<scream::Real, PACKN> Pack;

  static void test_index () {
    
  }

  static void test_range () {
    
  }

#define test_pack_gen_assign_op_all(op)
#define test_pack_gen_bin_op_all(op)
#define test_pack_gen_unary_unary_fn(op)
#define test_mask_gen_bin_op_all(op)

  static void run () {
    test_pack_gen_assign_op_all(=);
    test_pack_gen_assign_op_all(+=);
    test_pack_gen_assign_op_all(-=);
    test_pack_gen_assign_op_all(*=);
    test_pack_gen_assign_op_all(/=);

    test_pack_gen_bin_op_all(+);
    test_pack_gen_bin_op_all(-);
    test_pack_gen_bin_op_all(*);
    test_pack_gen_bin_op_all(/);

    test_pack_gen_unary_unary_fn(abs);
    test_pack_gen_unary_unary_fn(exp);
    test_pack_gen_unary_unary_fn(log);
    test_pack_gen_unary_unary_fn(log10);
    test_pack_gen_unary_unary_fn(tgamma);

    test_mask_gen_bin_op_all(==);
    test_mask_gen_bin_op_all(>=);
    test_mask_gen_bin_op_all(<=);
    test_mask_gen_bin_op_all(>);
    test_mask_gen_bin_op_all(<);

    test_index();
    test_range();
  }
};

TEST_CASE("Pack", "scream::pack") {
  TestPack<1>::run();
  TestPack<2>::run();
  TestPack<3>::run();
  TestPack<4>::run();
  TestPack<8>::run();
  TestPack<16>::run();
  TestPack<32>::run();
}

} // empty namespace
