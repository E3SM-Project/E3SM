#include <catch2/catch.hpp>
#include <numeric>

#include "share/expressions/base.hpp"
#include "share/expressions/binary_op.hpp"
#include "share/expressions/compare.hpp"
#include "share/expressions/conditional.hpp"
#include "share/expressions/field.hpp"

#include "share/expressions/evaluate.hpp"
#include "share/expressions/helpers.hpp"

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

TEST_CASE("expressions", "") {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::uniform_real_distribution<Real> pdf(0, 1);
  auto engine = setup_random_test();

  FieldLayout fl ({COL},{10});

  FieldIdentifier fid ("", fl, kg, "some_grid");
  Field f1 (fid.alias("f1"),true);
  Field f2 (fid.alias("f2"),true);
  randomize(f1, engine, pdf);
  randomize(f2, engine, pdf);

  SECTION ("field") {
    auto f1e = RealFieldExpression(f1);

    // Evaluate
    Field f3(fid.alias("f3"),true);
    f3.deep_copy(-1);
    evaluate(f1e,f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==v1h(i));
    }
  }
  SECTION ("f_plus_f") {
    auto sum = f1 + f2;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    evaluate(sum,f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==(v1h(i)+v2h(i)));
    }
  }
  SECTION ("f_plus_prod") {
    auto expression = f1 * f2 + f1;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    evaluate(expression,f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)*v2h(i) + v1h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
  SECTION ("prod_plus_div") {
    auto expression = f1 * f2 + f1 / f2;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    evaluate(expression,f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)*v2h(i) + v1h(i)/v2h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
  SECTION ("prod_plus_div_minus_f") {
    auto expression = (f1 * f2 + f1 / f2) - f1;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    evaluate(expression,f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)*v2h(i) + v1h(i)/v2h(i) - v1h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
  SECTION ("cmp") {
    auto f1e = RealFieldExpression(f1);
    auto f2e = RealFieldExpression(f2);
    auto cmp = f1e>=f2e;

    // Evaluate
    Field f3(fid.alias("f3"),true);
    evaluate(cmp,f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      REQUIRE (v3h(i)==static_cast<Real>(v1h(i)>=v2h(i)));
    }
  }
  SECTION ("conditional") {
    auto f1e = RealFieldExpression(f1);
    auto f2e = RealFieldExpression(f2);
    auto expression = conditional(f1e>=f2e,f1+f2,f1-f2);

    // Evaluate
    Field f3(fid.alias("f3"),true);
    evaluate(expression,f3);

    f3.sync_to_host();
    auto v1h = f1.get_view<const Real*,Host>();
    auto v2h = f2.get_view<const Real*,Host>();
    auto v3h = f3.get_view<const Real*,Host>();
    for (int i=0; i<fl.size(); ++i) {
      auto tgt = v1h(i)>=v2h(i) ? v1h(i)+v2h(i) : v1h(i)-v2h(i);
      REQUIRE (v3h(i)==tgt);
    }
  }
}

} // namespace scream
