#include <catch2/catch.hpp>
#include <numeric>

#include "share/field_math/expression.hpp"
#include "share/field_math/sum_expression.hpp"
#include "share/field_math/cmp_expression.hpp"
#include "share/field_math/field_evaluate.hpp"

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

TEST_CASE("expressions", "") {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using KT = Field::kt_dev;

  std::uniform_real_distribution<Real> pdf(0, 1);
  auto engine = setup_random_test();

  FieldLayout fl ({COL},{10});

  FieldIdentifier fid ("", fl, kg, "some_grid");
  Field f1 (fid.alias("f1"),true);
  Field f2 (fid.alias("f2"),true);
  randomize(f1, engine, pdf);
  randomize(f2, engine, pdf);

  // Setup the evaluation of f1+f2
  auto f1e = RealFieldEvaluate(f1);
  auto f2e = RealFieldEvaluate(f2);

  SECTION ("sum") {
    auto sum = f1e + f2e;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    auto v3 = f3.get_view<Real*>();
    auto lambda = KOKKOS_LAMBDA(int i) {
      v3(i) = sum(i);
    };
    auto policy = KT::RangePolicy(0,fl.size());
    Kokkos::parallel_for(policy,lambda);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==(v1h(i)+v2h(i)));
    }
  }
  SECTION ("cmp") {
    auto cmp = f1e>=f2e;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    auto v3 = f3.get_view<Real*>();
    auto lambda = KOKKOS_LAMBDA(int i) {
      v3(i) = cmp(i) ? f1e(i) : f2e(i);
    };
    auto policy = KT::RangePolicy(0,fl.size());
    Kokkos::parallel_for(policy,lambda);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==std::max(v1h(i),v2h(i)));
    }
  }
}

} // namespace scream
